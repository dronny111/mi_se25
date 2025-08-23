#!/usr/bin/env python3
"""
Classify CRISPR guides as inside/outside the bZIP‑1 domain of SlAREB1.

The script first tries to fetch the bZIP‑1 domain (InterPro IPR011616)
from the UniProt entry Q0PN11.  If UniProt does not list that domain
(e.g. because the annotation changed), you can supply the protein
coordinates manually with ``--domain-start`` and ``--domain-end``.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO

# -----------------------------------------------------------------
# 1)  UniProt helpers – now match only the InterPro ID (IPR011616)
# -----------------------------------------------------------------
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/"

def fetch_uniprot_entry(accession: str) -> dict:
    """Download the UniProt entry as JSON; abort with a clear message if it fails."""
    url = f"{UNIPROT_API}{accession}?format=json"
    r = requests.get(url)
    if r.status_code != 200:
        sys.exit(f"❌  Unable to fetch UniProt entry {accession} (HTTP {r.status_code})")
    return r.json()


def get_bzip1_coords(uniprot_json: dict) -> tuple[int, int] | None:
    """
    Return the (protein_start, protein_end) (1‑based inclusive) of the
    InterPro feature IPR011616 (bZIP‑1).  If the feature is not present,
    return ``None`` – the caller can then fall back to a user‑supplied
    interval.
    """
    for ft in uniprot_json.get("features", []):
        # We only need the InterPro identifier; description may vary.
        interpro = ft.get("interpro") or {}
        if interpro.get("id") == "IPR011616":
            # UniProt locations are 1‑based inclusive.
            start = int(ft["location"]["start"]["value"])
            end   = int(ft["location"]["end"]["value"])
            return start, end
    # Feature not found → signal the caller
    return None


# -----------------------------------------------------------------
# 2)  Convert protein aa positions → CDS nucleotide positions
# -----------------------------------------------------------------
def protein_to_cds_nucs(prot_start: int, prot_end: int) -> tuple[int, int]:
    """
    Convert protein aa positions (1‑based) to CDS nucleotide positions
    (1‑based, inclusive).  The conversion is simply:
        nuc_start = (prot_start‑1) * 3 + 1
        nuc_end   = prot_end * 3
    """
    nuc_start = (prot_start - 1) * 3 + 1
    nuc_end   = prot_end * 3
    return nuc_start, nuc_end


# -----------------------------------------------------------------
# 3)  Core routine – classify guides
# -----------------------------------------------------------------
def classify_guides(guides_csv: Path,
                    cds_fasta: Path,
                    uniprot_accession: str,
                    manual_prot_start: int | None,
                    manual_prot_end: int | None) -> None:
    """
    Reads the guide CSV produced by ``src/guide.py``, determines the
    bZIP‑1 domain (either from UniProt or from the user‑provided
    coordinates), and writes two CSV files:

        <guides>_in_domain.csv
        <guides>_out_of_domain.csv
    """
    # -----------------------------------------------------------------
    # 3.1  Load guide table
    # -----------------------------------------------------------------
    guides = pd.read_csv(guides_csv)

    # -----------------------------------------------------------------
    # 3.2  Load the CDS (just a sanity‑check on length)
    # -----------------------------------------------------------------
    cds_record = SeqIO.read(str(cds_fasta), "fasta")
    cds_len = len(cds_record.seq)
    print(f"🧬  CDS length = {cds_len:,} bp (from {cds_fasta.name})")

    # -----------------------------------------------------------------
    # 3.3  Determine the protein coordinates of the bZIP‑1 domain
    # -----------------------------------------------------------------
    prot_start, prot_end = None, None

    # 3.3.1  Try UniProt first
    uniprot_json = fetch_uniprot_entry(uniprot_accession)
    coords = get_bzip1_coords(uniprot_json)
    if coords:
        prot_start, prot_end = coords
        print(f"🔎  bZIP‑1 domain found in UniProt: protein aa {prot_start}-{prot_end}")
    else:
        print("⚠️  bZIP‑1 (IPR011616) NOT found in the UniProt entry.")
        # 3.3.2  Fall back to manual coordinates, if the user supplied them
        if manual_prot_start is not None and manual_prot_end is not None:
            prot_start, prot_end = manual_prot_start, manual_prot_end
            print(f"🔧  Using user‑supplied protein coordinates: aa {prot_start}-{prot_end}")
        else:
            sys.exit(
                "❌  No domain coordinates available.  Either:\n"
                "        • Update the UniProt entry (it may have been re‑annotated), or\n"
                "        • Supply the protein coordinates manually with\n"
                "          --domain-start <int> --domain-end <int>"
            )

    # -----------------------------------------------------------------
    # 3.4  Convert protein → nucleotide interval on the CDS
    # -----------------------------------------------------------------
    nuc_start, nuc_end = protein_to_cds_nucs(prot_start, prot_end)
    print(f"🧩  Corresponding CDS nucleotides: {nuc_start}-{nuc_end}")

    if nuc_end > cds_len:
        sys.exit(
            f"❌  Domain nucleotides ({nuc_end}) exceed the CDS length ({cds_len}).\n"
            "    Check that the supplied protein coordinates belong to this isoform."
        )

    # -----------------------------------------------------------------
    # 3.5  Overlap test (half‑open intervals, 0‑based)
    # -----------------------------------------------------------------
    def overlaps_domain(row):
        prot_start0 = row["genomic_start"]           # 0‑based start of the 20‑mer
        prot_end0   = row["genomic_end"] - 3        # drop the 3‑base PAM
        # domain as 0‑based half‑open interval:
        dom_start0 = nuc_start - 1
        dom_end0   = nuc_end
        return not (prot_end0 <= dom_start0 or prot_start0 >= dom_end0)

    guides["in_bzip1"] = guides.apply(overlaps_domain, axis=1)

    # -----------------------------------------------------------------
    # 3.6  Write the two output tables
    # -----------------------------------------------------------------
    base = guides_csv.with_suffix("")
    in_path  = Path(f"{base}_in_domain.csv")
    out_path = Path(f"{base}_out_of_domain.csv")

    guides[guides["in_bzip1"]].to_csv(in_path, index=False)
    guides[~guides["in_bzip1"]].to_csv(out_path, index=False)

    print(f"✅  {guides['in_bzip1'].sum()} guides **inside** the bZIP‑1 domain → {in_path}")
    print(f"✅  {len(guides) - guides['in_bzip1'].sum()} guides **outside** the domain → {out_path}")


# -----------------------------------------------------------------
# 4)  Command‑line interface
# -----------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Classify CRISPR guides as inside/outside the bZIP‑1 domain of SlAREB1."
    )
    parser.add_argument(
        "--guides", type=Path, required=True,
        help="CSV produced by src/guide.py (must contain genomic_start/end)."
    )
    parser.add_argument(
        "--cds", type=Path, required=True,
        help="FASTA of the SlAREB1 coding sequence (the same file used for guide generation)."
    )
    parser.add_argument(
        "--uniprot", type=str, default="Q0PN11",
        help="UniProt accession of the SlAREB1 protein (default Q0PN11)."
    )
    parser.add_argument(
        "--domain-start", type=int, default=None,
        help="Protein‑level start (aa) of the bZIP‑1 domain – used only if UniProt lookup fails."
    )
    parser.add_argument(
        "--domain-end", type=int, default=None,
        help="Protein‑level end (aa) of the bZIP‑1 domain – used only if UniProt lookup fails."
    )
    args = parser.parse_args()

    # Simple sanity check: if the user supplies one of the manual coordinates,
    # they must supply *both* so that the interval is well defined.
    if (args.domain_start is None) ^ (args.domain_end is None):
        sys.exit("❌  Either provide BOTH --domain-start AND --domain-end, or provide NEITHER.")

    classify_guides(
        guides_csv=args.guides,
        cds_fasta=args.cds,
        uniprot_accession=args.uniprot,
        manual_prot_start=args.domain_start,
        manual_prot_end=args.domain_end,
    )
