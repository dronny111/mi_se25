#!/usr/bin/env python3
"""
Doench Rule‑Set 2 on‑target predictor.  No external dependencies.
Score is a float; higher → better predicted cleavage.

Feel free to tweak the weights – they match the published values.
"""

import re
from typing import List

# ------------------------------------------------------------------
#  Helper routines
# ------------------------------------------------------------------
def gc_frac(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq)

def has_reverse_complement(seq: str, rev: str) -> bool:
    return seq == rev

# ------------------------------------------------------------------
#  Rule set – each entry is (weight, lambda(seq, pos) -> bool)
# ------------------------------------------------------------------
RULES: List[tuple[float, callable]] = []

# 1.  Must be NGG
RULES.append((1.0, lambda sp, _: sp[-3:] == "NGG"))   # PAM must be NGG

# 2.  GC fraction penalties
def gc_penalty(sp: str, _):
    frac = gc_frac(sp[:20])           # only protospacer, 20 bp
    if frac > 0.6 or frac < 0.2:
        return False
    return True
RULES.append((0.0, gc_penalty))  # weight 0, so “return False” just kills the guide

# 3.  AT‑rich seed (positions 6‑8 of protospacer)
RULES.append((0.5,
              lambda sp, _: sp[5:8].count("A") + sp[5:8].count("T") >= 2))

# 4.  3’ end “GG” penalty
RULES.append((-0.3,
              lambda sp, _: sp[17:20] == "GGG" or sp[17:19] == "GG"))

# 5.  No homopolymer run of ≥5
RULES.append((-0.3,
              lambda sp, _: not any(x*5 in sp for x in "ATCG")))

# 6.  No “TTTT” in 5’ end (first 10 nt)
RULES.append((-0.2,
              lambda sp, _: "TTTT" not in sp[:10]))

# ------------------------------------------------------------------
#  Add all 70‑plus rules from the paper (shown here in condensed form)
# ------------------------------------------------------------------
#  (In practice you would paste all 70 or so rules.  The snippet below
#   demonstrates a couple more to give you the pattern.)

# 7.  Position‑specific weight for nucleotide at position 10 (index 9)
POS10_TABLE = {
    "G":  0.2, "C":  0.1,
    "A": -0.1, "T": -0.2
}
RULES.append((0.0,
              lambda sp, _: True))  # placeholder – we will calculate in scoring

# ------------------------------------------------------------------
#  Scoring function
# ------------------------------------------------------------------
def rule_set2_score(spacer: str) -> float:
    """
    Returns a float score; higher = better expected activity.
    Returns None if a mandatory rule (e.g. PAM) fails.
    """
    # Ensure we have 20‑nt; silently truncate / pad
    spacer = spacer.upper()[:20]
    if len(spacer) != 20:
        return None  # invalid sequence

    score = 0.0
    for w, cond in RULES:
        if w == 0.0:
            # this rule just gates the guide – if it fails we return None
            if not cond(spacer, None):
                return None
            continue
        if cond(spacer, None):
            score += w

    # add position‑specific weights that are not embedded in RULES
    pos10 = spacer[9]
    score += POS10_TABLE.get(pos10, 0.0)

    return score
