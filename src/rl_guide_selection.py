#!/usr/bin/env python3
import json, random
import numpy as np
import pandas as pd
from typing import List, Tuple

def parse_outcome_json(s: str) -> List[Tuple[int,float]]:
    try:
        lst = json.loads(s)
    except:
        return [(+1,1.0)]
    return [(item["size"], item["prob"]) for item in lst]

def simulate_reward(protospacer, rule2_score, outcome_json) -> float:
    outcomes = parse_outcome_json(outcome_json)
    sizes, probs = zip(*outcomes)
    size = random.choices(sizes, weights=probs, k=1)[0]
    frameshift = (size % 3) != 0
    return float(frameshift) * rule2_score

class ThompsonBandit:
    def __init__(self, n):
        self.n = n
        self.alpha = np.ones(n)
        self.beta  = np.ones(n)

    def select_arm(self):
        return int(np.argmax(np.random.beta(self.alpha, self.beta)))

    def update(self, arm, reward):
        if reward > 0.0:
            self.alpha[arm] += 1
        else:
            self.beta[arm]  += 1

def run_thompson_sampling(guides_df: pd.DataFrame,
                          n_pulls: int = 200,
                          top_n: int = 5,
                          random_seed: int = 42) -> pd.DataFrame:
    random.seed(random_seed); np.random.seed(random_seed)
    bandit = ThompsonBandit(len(guides_df))

    for _ in range(n_pulls):
        arm = bandit.select_arm()
        row = guides_df.iloc[arm]
        r = simulate_reward(row["protospacer"], row["score"], row["predicted_outcome"])
        bandit.update(arm, r)

    post = bandit.alpha / (bandit.alpha + bandit.beta)
    guides_df = guides_df.copy()
    guides_df["alpha"] = bandit.alpha
    guides_df["beta"]  = bandit.beta
    guides_df["posterior_mean"] = post
    top = guides_df.nlargest(top_n, "posterior_mean").reset_index(drop=True)
    return top
