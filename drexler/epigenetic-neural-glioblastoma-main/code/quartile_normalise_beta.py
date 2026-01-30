#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import sys


def read_single_sample_beta(path: str) -> pd.Series:
    # First line: sample id
    with open(path, "r") as f:
        first = f.readline().strip()

    if not first:
        raise ValueError("Input file is empty or missing sample ID on first line.")

    sample_id = first

    # Remaining lines: cpg,beta
    df = pd.read_csv(path, skiprows=1, header=None, names=["cpg", "beta"])
    if df.shape[1] != 2:
        raise ValueError("Expected two columns (cpg,beta) after the first line.")

    df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
    s = df.set_index("cpg")["beta"]
    s.name = sample_id
    return s


def beta_to_m(beta: pd.Series, eps=1e-6) -> pd.Series:
    beta = beta.clip(eps, 1 - eps)
    return np.log2(beta / (1 - beta))


def m_to_beta(M: pd.Series) -> pd.Series:
    x = np.power(2.0, M)
    return x / (1 + x)


def reference_quantile_map_single_sample(m_sample: pd.Series, m_ref: pd.DataFrame) -> pd.Series:
    """
    Map sample M-values to the reference distribution using rank-based mapping.
    m_ref: rows=CpGs, cols=reference samples (Illumina/training), M-values
    """
    # Use only CpGs shared between sample and reference
    common = m_sample.index.intersection(m_ref.index)
    if len(common) < 1000:
        raise ValueError(f"Too few overlapping CpGs with reference:

