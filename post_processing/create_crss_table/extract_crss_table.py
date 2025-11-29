
# standard lib
import os
import re
import json
import shutil
import tarfile
import tempfile

# 3rd party lib
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict


def main():
    with open("config.json", "r") as f:
        config = json.load(f)
    simulationDir = Path(config["data_path"])
    output_dir = Path(config["output_path"])

    df_data = Path("/scratch/hyunsol/MoDELib2_20251027_Fork_git_rotate/MoDELib2/post_processing/create_crss_table/output/partial_orient_615AA_crss_0-1_almg5_total/all_rate_data.csv")
    df = pd.read_csv(df_data)

    # generate CRSS table
    seeds = df["seed"].unique()
    threshold = 1e-10  # set your cut-off

    data_frames = []
    data_columns = [
        "d_type",
        "noise_type",
        "seed",
        "d_char",
        "alloy",
        "length",
        "str_step_size",
        "crss",
    ]
    for seed in seeds:
        local_df = df[df["seed"] == seed]
        high_df = local_df[local_df["dydx_mean"] > threshold]  # 1. apply threshold
        crss_above_row = high_df.loc[high_df["dydx_mean"].idxmin()]
        crss_above = crss_above_row["stress"]

        # rows whose stress is below the CRSS
        below = local_df[local_df["stress"] < crss_above]
        if not len(below):
            crss = crss_above
            crss_search_step = crss_above
        else:
            # pick the one with the highest such stress
            entry_below = below.loc[below["stress"].idxmax(), "stress"]
            crss_search_step = crss_above - entry_below
            crss = crss_above - (crss_search_step / 2)
        new_entry = [
            (
                crss_above_row["d_type"],
                crss_above_row["noise_type"],
                crss_above_row["seed"],
                crss_above_row["d_char"],
                crss_above_row["alloy"],
                crss_above_row["length"],
                crss_search_step,
                crss,
            )
        ]
        data_frames.append(pd.DataFrame(new_entry, columns=data_columns))

    df_new = pd.concat(data_frames, ignore_index=True)
    df_new = df_new.sort_values(by=["seed"])
    df_new.to_csv(output_dir / "all_crss_data.csv", index=False)

if __name__ == "__main__":
    main()
