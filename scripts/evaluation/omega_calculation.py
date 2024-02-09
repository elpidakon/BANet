# Kontsioti, Maskell, Anderson & Pirmohamed, Identifying drug-drug interactions in spontaneous reports utilizing signal detection and biological plausibility aspects (2024)
# Script for calculation of Omega values for the positive and negative controls included in the CRESCENDDI reference set.

import pandas as pd
import math
import numpy as np
from scipy import stats
import os

def omega_calculation(datafile, output):
    df = pd.read_excel(datafile, header=0)
    df["omega_025"] = ""
    for index, row in df.iterrows():
        if 0 not in {
            row["d1_d2_counter"],
            row["not_d1_d2_counter"],
            row["d1_not_d2_counter"],
            row["n_111"],
        }:
            f_00 = row["n_001"] / row["not_d1_not_d2_counter"]
            f_01 = row["n_011"] / row["not_d1_d2_counter"]
            f_10 = row["n_101"] / row["d1_not_d2_counter"]
            f_11 = row["n_111"] / row["d1_d2_counter"]
            g_11 = 1 - (
                1
                / (
                    max(f_00 / (1 - f_00), f_10 / (1 - f_10))
                    + max(f_00 / (1 - f_00), f_01 / (1 - f_01))
                    - f_00 / (1 - f_00)
                    + 1
                )
            )
            E_111 = g_11 * row["d1_d2_counter"]
            mu = stats.gamma.ppf(0.025, row["n_111"] + 0.5, scale=1 / (E_111 + 0.5))
            omega_025 = math.log2(mu)
            df.loc[index, 'omega_025'] = omega_025
        else:
            df.loc[index, "omega_025"] = np.nan
    df.to_excel(output, index=False)


omega_calculation(
    "data/evaluation/DR1_FAERS_COUNTS.xlsx",
    "output/evaluation/DR1_OMEGA_VALUES.xlsx",
)
omega_calculation(
    "data/evaluation/DR2_FAERS_COUNTS.xlsx",
    "output/evaluation/DR2_OMEGA_VALUES.xlsx",
)
