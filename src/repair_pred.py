#!/usr/bin/env python3
"""
Micro‑homology‑mediated deletion (MHD) repair outcome predictor.

Given a 20‑nt protospacer (the sequence that the sgRNA pairs with) we:

1. Assume the Cas9 cut occurs **3 bp upstream of the PAM** on the target strand,
   i.e. between bases 17 and 18 of the protospacer (0‑based indexing).

2. Scan a 30‑bp window centred on the cut site (15 bp upstream, 15 bp downstream)
   for all *identical* k‑mers (k ≥ 2).  Each pair of matching k‑mers defines a
   potential micro‑homology tract that can be used by the cell to ligate the
   broken ends, producing a deletion whose length is the distance between the
   two repeats **plus** the length of the repeat itself.

3. Assign a simple “strength” to each tract = *k* (the longer the repeat,
   the more likely it is to be used).  The raw strength is turned into a
   probability by normalising over all candidate tracts.

4. Add a **baseline 10 % probability for a +1 insertion** (the most common
   non‑deletion outcome observed for SpCas9).

The function returns a **list of (description, size, probability)** tuples,
sorted from most to least likely.  Example output for a guide:

    [('del‑8 (micro‑homology 4 bp)', -8, 0.46),
     ('del‑3 (micro‑homology 3 bp)', -3, 0.34),
     ('+1 insertion', +1, 0.10)]

The description is purposely concise because the list is written into the CSV
produced by ``src/guide.py`` (column ``predicted_outcome``).  Downstream
analyses can pick the top entry or keep the full list as a JSON string.
"""

from __future__ import annotations
import json
from typing import List, Tuple


CUT_POS_IN_PROTOSPACER = 17          # 0‑based; cut is between 17 & 18


def _microhomology_candidates(protospacer: str) -> List[Tuple[int, int, int]]:
    """
    Find all micro‑homology tracts (k ≥ 2) that straddle the cut site.

    Returns a list of tuples ``(repeat_len, del_len, start_upstream)`` where:
        * ``repeat_len`` – length of the identical repeat (k)
        * ``del_len``    – total number of bases removed (negative number)
        * ``start_upstream`` – index (0‑based) of the upstream copy of the repeat
    """
    # Build a 30‑bp window centred on the cut site (if the protospacer is
    # shorter at either side we simply use the available bases).
    left  = protospacer[:CUT_POS_IN_PROTOSPACER]          # bases 0 … 16
    right = protospacer[CUT_POS_IN_PROTOSPACER:]          # bases 17 … 19

    # We will compare every possible k‑mer in the left side with every k‑mer
    # in the right side.  Only exact matches are considered.
    candidates = []
    max_k = min(len(left), len(right))   # cannot be larger than the shorter side
    for k in range(2, max_k + 1):        # start at 2‑bp repeat
        # slide a window of length k over the left side
        for i in range(0, len(left) - k + 1):
            left_kmer = left[i:i + k]
            # slide over the right side and look for the same k‑mer
            for j in range(0, len(right) - k + 1):
                if left_kmer == right[j:j + k]:
                    # Deletion length = (distance between the two copies) + k
                    # distance = (len(left)-i) + j   (bases that will be removed)
                    del_len = - ( (len(left) - i) + j + k )
                    candidates.append( (k, del_len, i) )
    return candidates


def predict_repair_outcome(protospacer: str) -> List[Tuple[str, int, float]]:
    """
    Main entry point – given a 20‑nt protospacer returns a ranked list of
    predicted repair outcomes.

    Returns
    -------
    List[ ( description:str , size:int , probability:float ) ]
        * ``size`` is **negative** for deletions, **positive** for insertions.
        * probabilities sum to 1.0.
    """
    # -------------------------------------------------
    # 1)  Gather micro‑homology candidates
    # -------------------------------------------------
    candidates = _microhomology_candidates(protospacer)

    # -------------------------------------------------
    # 2)  Convert raw micro‑homology “strength” into probabilities
    # -------------------------------------------------
    # Strength = repeat length (k); we keep it linear.
    if candidates:
        total_strength = sum(k for k, _, _ in candidates)
        mh_probs = [
            (k, del_len, k / total_strength)   # (repeat_len, del_len, prob)
            for k, del_len, _ in candidates
        ]
    else:
        mh_probs = []   # no micro‑homology tracts found

    # -------------------------------------------------
    # 3)  Add a baseline +1 insertion (10 % of the total probability mass)
    # -------------------------------------------------
    INSERTION_PROB = 0.10
    # Reduce the micro‑homology probabilities proportionally so that total = 0.90
    scaling = 1.0 - INSERTION_PROB
    mh_probs = [(k, d, p * scaling) for k, d, p in mh_probs]

    # -------------------------------------------------
    # 4)  Build the final list of (desc, size, prob)
    # -------------------------------------------------
    outcome_list: List[Tuple[str, int, float]] = []

    for repeat_len, del_len, prob in mh_probs:
        desc = f"del{del_len} (micro‑homology {repeat_len} bp)"
        outcome_list.append( (desc, del_len, prob) )

    # Insert the +1 insertion entry (if we have any micro‑homology entries,
    # its probability will be the 10 % we reserved above; otherwise it is 1.0)
    insertion_prob = INSERTION_PROB if mh_probs else 1.0
    outcome_list.append( ("+1 insertion", +1, insertion_prob) )

    # -------------------------------------------------
    # 5)  Normalise (guard against rounding errors) and sort
    # -------------------------------------------------
    total = sum(p for _, _, p in outcome_list)
    outcome_list = [(d, s, p / total) for d, s, p in outcome_list]
    outcome_list.sort(key=lambda x: x[2], reverse=True)   # most likely first

    return outcome_list


def serialize_outcome(outcome: List[Tuple[str, int, float]]) -> str:
    """
    Helper that turns the list of tuples into a compact JSON string – useful
    for writing to a CSV column.
    """
    # Convert to a list of dicts for a more readable JSON payload
    dlist = [{"desc": d, "size": s, "prob": round(p, 4)} for d, s, p in outcome]
    return json.dumps(dlist, separators=(",", ":"))
