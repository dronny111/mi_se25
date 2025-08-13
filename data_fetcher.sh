#!/usr/bin/env bash
set -euo pipefail

# Usage:
# 1) Save as fetch_SlAREB1_SGN_first.sh
# 2) chmod +x fetch_SlAREB1_SGN_first.sh
# 3) ./fetch_SlAREB1_SGN_first.sh
#
# Requirements:
# - curl
# - jq
# - samtools (optional but recommended for promoter extraction if using local FASTA)
# - zip
# - efetch/edirect (optional) or curl to NCBI efetch
#
# This script tries SGN first for the versioned Solyc ID (Solyc04g078840.*).
# If SGN direct fasta links are unavailable, it falls back to Ensembl REST.
# Also downloads AY530758 from NCBI.
#
# If any step fails due to server-side HTML changes on SGN, the script will
# report the step and create the manifest entries it could.

OUTDIR="SlAREB1_bundle_SGN_first"
ZIPNAME="SlAREB1_sequences_bundle_SGN_first.zip"
ENSEMBL_REST="https://rest.ensembl.org"
NCBI_EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
GENE_BASE="Solyc04g078840"
SPECIES="solanum_lycopersicum"

mkdir -p "${OUTDIR}"/{genomic,cDNA,CDS,protein,promoter,homologs}
manifest="${OUTDIR}/manifest.tsv"
readme="${OUTDIR}/README.txt"

echo -e "file_name\tcultivar/species\tsequence_type\tlength_bp\taccession_or_id\tcoords\tdatabase\tsource_version\tassembly_year\turl" > "$manifest"

cat > "$readme" <<EOF
SlAREB1 bundle (SGN-first) README
---------------------------------
This bundle attempts to download SlAREB1 (Solyc04g078840.*) sequences from SGN first.
If direct SGN fasta downloads are not available the script will fallback to Ensembl REST.

Folders:
- genomic/, cDNA/, CDS/, protein/, promoter/, homologs/
Naming:
- SlAREB1_<CultivarOrSpecies>_<sequenceType>.fasta

Manifest contains direct URLs where possible.
EOF

# helper to check HTTP 200
http_ok() {
  local url=$1
  status=$(curl -s -o /dev/null -w "%{http_code}" "$url")
  [ "$status" = "200" ]
}

# 1) Try SGN for versioned IDs .4 .3 .2 .1 (descending -> prefer higher)
SGN_BASE="https://solgenomics.net/locus"
found_sgn_id=""
for v in 5 4 3 2 1; do
  candidate="${GENE_BASE}.${v}"
  url="${SGN_BASE}/${candidate}/view"
  echo "Checking SGN for ${candidate}..."
  if http_ok "$url"; then
    echo "Found SGN page: $url"
    found_sgn_id="$candidate"
    found_sgn_url="$url"
    break
  fi
done

if [ -z "$found_sgn_id" ]; then
  echo "No SGN page found for ${GENE_BASE} with typical suffixes. Will fallback to Ensembl symbol searches."
else
  echo "Using SGN gene id: $found_sgn_id"
fi

# 2) If SGN page found, attempt to locate direct FASTA links on SGN page
#    (SGN pages often include a "Download sequence" link. We'll try to scrape a href with "download" or ".fa" or "fasta").
downloaded_from_sgn=false
if [ -n "${found_sgn_id:-}" ]; then
  echo "Attempting to scrape SGN page for download links..."
  page_html=$(curl -s "$found_sgn_url")
  # Look for direct links ending in .fa .fasta .txt or containing "download"
  # This is best-effort; SGN page structure may vary.
  links=$(echo "$page_html" | grep -Eo 'href="[^"]+"' | sed 's/href="//;s/"$//' | grep -E '\.fa(sta)?$|download|/download' || true)
  if [ -n "$links" ]; then
    echo "Found candidate download links on SGN page. Trying them..."
    for l in $links; do
      # make URL absolute if needed
      if [[ "$l" =~ ^/ ]]; then
        l="https://solgenomics.net${l}"
      fi
      echo "Trying $l"
      if http_ok "$l"; then
        # try to infer type by filename
        fname=$(basename "$l")
        # handle common names
        if echo "$fname" | grep -Ei "genome|genomic|dna" >/dev/null 2>&1; then
          out="${OUTDIR}/genomic/SlAREB1_SGN_${found_sgn_id}_genomic.fasta"
        elif echo "$fname" | grep -Ei "cdna|mrna|transcript" >/dev/null 2>&1; then
          out="${OUTDIR}/cDNA/SlAREB1_SGN_${found_sgn_id}_cDNA.fasta"
        elif echo "$fname" | grep -Ei "cds" >/dev/null 2>&1; then
          out="${OUTDIR}/CDS/SlAREB1_SGN_${found_sgn_id}_CDS.fasta"
        elif echo "$fname" | grep -Ei "pep|prot|protein" >/dev/null 2>&1; then
          out="${OUTDIR}/protein/SlAREB1_SGN_${found_sgn_id}_protein.fasta"
        else
          out="${OUTDIR}/genomic/SlAREB1_SGN_${found_sgn_id}_misc_$(basename "$l")"
        fi
        curl -s -L "$l" -o "$out"
        if [ -s "$out" ]; then
          llen=$(grep -v '^>' "$out" | tr -d '\n' | wc -c || echo 0)
          echo -e "$(basename "$out")\tSGN\tunknown\t${llen}\t${found_sgn_id}\t-\tSGN\tunknown\tNA\t${l}" >> "$manifest"
          downloaded_from_sgn=true
        fi
      fi
    done
  else
    echo "No direct fasta/download links detected on SGN page (or page structure incompatible). Will fallback to Ensembl."
  fi
fi

# 3) If SGN download didn't produce the canonical sequences, fallback to Ensembl
if [ "$downloaded_from_sgn" = false ]; then
  echo "Falling back to Ensembl REST to fetch sequences using the symbol or versioned IDs."

  # Try Ensembl lookup with versioned ID first (if we found one)
  ensembl_gene_id=""
  if [ -n "${found_sgn_id:-}" ]; then
    echo "Trying Ensembl lookup with symbol ${found_sgn_id}..."
    gene_json=$(curl -s "${ENSEMBL_REST}/lookup/symbol/${SPECIES}/${found_sgn_id}?content-type=application/json")
    if ! echo "$gene_json" | jq -e '.error' >/dev/null 2>&1; then
      ensembl_gene_id=$(echo "$gene_json" | jq -r '.id')
    fi
  fi

  # If not found, try lookup by base symbol (no suffix)
  if [ -z "${ensembl_gene_id:-}" ]; then
    echo "Trying Ensembl lookup with symbol ${GENE_BASE}..."
    gene_json=$(curl -s "${ENSEMBL_REST}/lookup/symbol/${SPECIES}/${GENE_BASE}?content-type=application/json")
    if ! echo "$gene_json" | jq -e '.error' >/dev/null 2>&1; then
      ensembl_gene_id=$(echo "$gene_json" | jq -r '.id')
    fi
  fi

  if [ -z "${ensembl_gene_id:-}" ]; then
    echo "Ensembl lookup failed for both ${found_sgn_id:-<none>} and ${GENE_BASE}. At this point you can:"
    echo " - provide the exact versioned Solyc ID (e.g., Solyc04g078840.2), OR"
    echo " - download the sequences manually from SGN and place them in the respective folders under ${OUTDIR}/"
    echo "The script will still proceed to fetch AY530758 and orthologues via Ensembl where possible."
  else
    echo "Ensembl gene id found: $ensembl_gene_id"
    # fetch gene JSON to get coords and assembly
    chr=$(echo "$gene_json" | jq -r '.seq_region_name // "NA"')
    start=$(echo "$gene_json" | jq -r '.start // "NA"')
    end=$(echo "$gene_json" | jq -r '.end // "NA"')
    strand=$(echo "$gene_json" | jq -r '.strand // 1')
    assembly_name=$(echo "$gene_json" | jq -r '.assembly_name // "unknown"')

    # fetch genomic sequence
    echo "Fetching genomic sequence..."
    curl -s "${ENSEMBL_REST}/sequence/id/${ensembl_gene_id}?type=genomic;content-type=text/plain" -o "${OUTDIR}/genomic/SlAREB1_Ensembl_${ensembl_gene_id}_genomic.fasta"
    glen=$(grep -v '^>' "${OUTDIR}/genomic/SlAREB1_Ensembl_${ensembl_gene_id}_genomic.fasta" | tr -d '\n' | wc -c || echo 0)
    echo -e "SlAREB1_Ensembl_${ensembl_gene_id}_genomic.fasta\tHeinz1706\tgenomic\t${glen}\t${ensembl_gene_id}\tchr${chr}:${start}-${end}\tEnsembl\t${assembly_name}\tNA\t${ENSEMBL_REST}/sequence/id/${ensembl_gene_id}?type=genomic" >> "$manifest"

    # fetch transcripts (cdna, cds, protein)
    trans_json=$(curl -s "${ENSEMBL_REST}/overlap/id/${ensembl_gene_id}?feature=transcript;content-type=application/json")
    trans_ids=$(echo "$trans_json" | jq -r '.[].id' | sort -u || true)
    for tid in $trans_ids; do
      curl -s "${ENSEMBL_REST}/sequence/id/${tid}?type=cdna;content-type=text/plain" -o "${OUTDIR}/cDNA/SlAREB1_Ensembl_${tid}_cDNA.fasta"
      curl -s "${ENSEMBL_REST}/sequence/id/${tid}?type=cds;content-type=text/plain" -o "${OUTDIR}/CDS/SlAREB1_Ensembl_${tid}_CDS.fasta"
      curl -s "${ENSEMBL_REST}/sequence/id/${tid}?type=protein;content-type=text/plain" -o "${OUTDIR}/protein/SlAREB1_Ensembl_${tid}_protein.fasta"
      l1=$(grep -v '^>' "${OUTDIR}/cDNA/SlAREB1_Ensembl_${tid}_cDNA.fasta" | tr -d '\n' | wc -c || echo 0)
      l2=$(grep -v '^>' "${OUTDIR}/CDS/SlAREB1_Ensembl_${tid}_CDS.fasta" | tr -d '\n' | wc -c || echo 0)
      l3=$(grep -v '^>' "${OUTDIR}/protein/SlAREB1_Ensembl_${tid}_protein.fasta" | tr -d '\n' | wc -c || echo 0)
      echo -e "SlAREB1_Ensembl_${tid}_cDNA.fasta\tHeinz1706\tcDNA\t${l1}\t${tid}\t-\tEnsembl\t${assembly_name}\tNA\t${ENSEMBL_REST}/sequence/id/${tid}?type=cdna" >> "$manifest"
      echo -e "SlAREB1_Ensembl_${tid}_CDS.fasta\tHeinz1706\tCDS\t${l2}\t${tid}\t-\tEnsembl\t${assembly_name}\tNA\t${ENSEMBL_REST}/sequence/id/${tid}?type=cds" >> "$manifest"
      echo -e "SlAREB1_Ensembl_${tid}_protein.fasta\tHeinz1706\tprotein\t${l3}\t${tid}\t-\tEnsembl\t${assembly_name}\tNA\t${ENSEMBL_REST}/sequence/id/${tid}?type=protein" >> "$manifest"
    done

    # promoter extraction via Ensembl region sequence endpoint (2kb upstream)
    PROM_LEN=2000
    if [ "$strand" -eq 1 ]; then
      prom_start=$((start - PROM_LEN))
      if [ $prom_start -lt 1 ]; then prom_start=1; fi
      prom_end=$((start - 1))
    else
      prom_start=$((end + 1))
      prom_end=$((end + PROM_LEN))
    fi
    prom_url="${ENSEMBL_REST}/sequence/region/${SPECIES}/${chr}:${prom_start}..${prom_end}:1?content-type=text/plain"
    echo "Fetching promoter via Ensembl: ${chr}:${prom_start}-${prom_end}..."
    curl -s -H "Content-Type: text/plain" "$prom_url" -o "${OUTDIR}/promoter/SlAREB1_Ensembl_promoter2kb.fasta"
    plen=$(grep -v '^>' "${OUTDIR}/promoter/SlAREB1_Ensembl_promoter2kb.fasta" | tr -d '\n' | wc -c || echo 0)
    echo -e "SlAREB1_Ensembl_promoter2kb.fasta\tHeinz1706\tpromoter\t${plen}\t-\tchr${chr}:${prom_start}-${prom_end}\tEnsembl\t${assembly_name}\tNA\t${prom_url}" >> "$manifest"
  fi
fi

# 4) Always fetch AY530758 (NCBI GenBank cDNA)
AYID="AY530758"
echo "Fetching AY530758 from NCBI..."
if command -v efetch >/dev/null 2>&1; then
  efetch -db nucleotide -id "$AYID" -format fasta > "${OUTDIR}/cDNA/SlAREB1_${AYID}_Moneymaker_cDNA.fasta"
else
  curl -s "${NCBI_EFETCH}?db=nucleotide&id=${AYID}&rettype=fasta&retmode=text" -o "${OUTDIR}/cDNA/SlAREB1_${AYID}_Moneymaker_cDNA.fasta"
fi
aylen=$(grep -v '^>' "${OUTDIR}/cDNA/SlAREB1_${AYID}_Moneymaker_cDNA.fasta" | tr -d '\n' | wc -c || echo 0)
echo -e "SlAREB1_${AYID}_Moneymaker_cDNA.fasta\tMoneymaker\tcDNA\t${aylen}\t${AYID}\t-\tNCBI\tnucleotide\tNA\thttps://www.ncbi.nlm.nih.gov/nuccore/${AYID}" >> "$manifest"

# 5) Fetch orthologues via Ensembl Compara (best-effort)
echo "Attempting to fetch orthologues via Ensembl Compara (if Ensembl gene ID known)..."
if [ -n "${ensembl_gene_id:-}" ]; then
  hom_url="${ENSEMBL_REST}/homology/id/${ensembl_gene_id}?content-type=application/json;type=orthologues;format=condensed"
  hom_json=$(curl -s "$hom_url")
  # attempt to download orthologues for a list of common Solanaceae species
  targets=("solanum_tuberosum" "capsicum_annuum" "solanum_melongena" "nicotiana_tabacum" "solanum_pimpinellifolium" "solanum_pennellii")
  for sp in "${targets[@]}"; do
    echo "Looking for orthologue in ${sp}..."
    entries=$(echo "$hom_json" | jq -r --arg sp "$sp" '.data[].homologies[] | select(.target.species==$sp) | .target.id' || true)
    if [ -z "$entries" ]; then
      echo "  No orthologue found in ${sp} via Ensembl."
      continue
    fi
    for tid in $entries; do
      echo "  Found orthologue ${tid} in ${sp} — fetching sequences..."
      curl -s "${ENSEMBL_REST}/sequence/id/${tid}?type=cdna;content-type=text/plain" -o "${OUTDIR}/homologs/SlAREB1_${sp}_transcript_${tid}_cDNA.fasta"
      curl -s "${ENSEMBL_REST}/sequence/id/${tid}?type=cds;content-type=text/plain" -o "${OUTDIR}/homologs/SlAREB1_${sp}_transcript_${tid}_CDS.fasta"
      curl -s "${ENSEMBL_REST}/sequence/id/${tid}?type=protein;content-type=text/plain" -o "${OUTDIR}/homologs/SlAREB1_${sp}_transcript_${tid}_protein.fasta"
      c1=$(grep -v '^>' "${OUTDIR}/homologs/SlAREB1_${sp}_transcript_${tid}_cDNA.fasta" | tr -d '\n' | wc -c || echo 0)
      c2=$(grep -v '^>' "${OUTDIR}/homologs/SlAREB1_${sp}_transcript_${tid}_CDS.fasta" | tr -d '\n' | wc -c || echo 0)
      c3=$(grep -v '^>' "${OUTDIR}/homologs/SlAREB1_${sp}_transcript_${tid}_protein.fasta" | tr -d '\n' | wc -c || echo 0)
      echo -e "SlAREB1_${sp}_transcript_${tid}_cDNA.fasta\t${sp}\tcDNA\t${c1}\t${tid}\t-\tEnsembl\tunknown\tNA\t${ENSEMBL_REST}/sequence/id/${tid}?type=cdna" >> "$manifest"
      echo -e "SlAREB1_${sp}_transcript_${tid}_CDS.fasta\t${sp}\tCDS\t${c2}\t${tid}\t-\tEnsembl\tunknown\tNA\t${ENSEMBL_REST}/sequence/id/${tid}?type=cds" >> "$manifest"
      echo -e "SlAREB1_${sp}_transcript_${tid}_protein.fasta\t${sp}\tprotein\t${c3}\t${tid}\t-\tEnsembl\tunknown\tNA\t${ENSEMBL_REST}/sequence/id/${tid}?type=protein" >> "$manifest"
    done
  done
else
  echo "Ensembl gene id unknown — skipping Ensembl orthologue fetch. You can re-run after providing versioned Solyc ID or after manually downloading SGN files."
fi

# 6) Finalize zip
echo "Creating zip archive ${ZIPNAME}..."
zip -r "$ZIPNAME" "$OUTDIR" >/dev/null
echo "Done. Archive: ${ZIPNAME}"
echo "Manifest: ${manifest}"
echo "If some SGN-specific files are missing, download them manually from the SGN locus page and place them under the appropriate folder in ${OUTDIR}/ and re-run 'zip -r' to refresh the archive."
