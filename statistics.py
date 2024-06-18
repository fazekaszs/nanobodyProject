import pandas as pd


def get_global_statistics(table: pd.DataFrame):

    table_stats = pd.concat({
        "Minimum": (t_min := table.min()),
        "Maximum": (t_max := table.max()),
        "Range": t_max - t_min,
        "Mean": table.mean(),
        "StD": table.std(),
        "Median": (t_med := table.median()),
        "MAD": (table - t_med).abs().median(),
        "Sarle\'s bimodality coefficient": (table.skew() ** 2 + 1) / (table.kurt() + 3),
        "Percentage greater than 0": table.map(lambda x: 100 * int(x > 0)).sum() / table.shape[0]
    }, axis=1)

    print(table_stats.to_string(), "\n")

    ordered_table_col_names = [
        "SWT:mCSM-PPI2", "SD1:mCSM-PPI2", "SD2:mCSM-PPI2",
        "SWT:MutaBind2", "SD1:MutaBind2", "SD2:MutaBind2",
        "SWT:SAAMBE-3D", "SD1:SAAMBE-3D", "SD2:SAAMBE-3D",
    ]
    ordered_table_row_names = ["H11H4", "Nb20", "ab8"]
    muts_best = list()
    muts_worst = list()

    for nanobody in ordered_table_row_names:

        muts_best_row = list()
        muts_worst_row = list()
        for col_name in ordered_table_col_names:  # S-RBD variants + predictors

            srbd_variant, predictor = col_name.split(":")
            relevant_values = table[predictor].filter(like=f"{nanobody}_{srbd_variant}:").sort_values()

            current_bests = map(lambda x: x.split(":")[1], relevant_values.index[:3])
            current_bests = ",".join(current_bests)
            muts_best_row.append(current_bests)

            current_worsts = map(lambda x: x.split(":")[1], relevant_values.index[-3:])
            current_worsts = ",".join(current_worsts)
            muts_worst_row.append(current_worsts)

        muts_best.append(muts_best_row)
        muts_worst.append(muts_worst_row)

    muts_best = pd.DataFrame(data=muts_best, index=ordered_table_row_names, columns=ordered_table_col_names)
    muts_worst = pd.DataFrame(data=muts_worst, index=ordered_table_row_names, columns=ordered_table_col_names)

    print(muts_best.to_csv(sep="\t"))
    print(muts_worst.to_csv(sep="\t"))


def get_z_statistics(table: pd.DataFrame):

    t_center = table.median()  # table.mean() or table.median()
    t_scale = (table - t_center).abs().median()  # table.std() or (table - t_center).abs().median()
    t_z = (table - t_center) / t_scale
    t_z = t_z.mean(axis=1)
    t_z = t_z.sort_values()

    nanobodies, variants = ["H11H4", "Nb20", "ab8"], ["SWT", "SD1", "SD2"]

    rows = list()
    values = list()

    for var in variants:
        for nb in nanobodies:

            key = f"{nb}_{var}"
            rows.append(key)

            current_muts = t_z.filter(like=key)
            print(f"Number of interface mutations for {key} is: {len(current_muts)}")
            current_muts = current_muts[:5]
            current_muts = [mut.replace(f"{key}:", "") for mut in current_muts.index]
            values.append(current_muts)

    mut_table = pd.DataFrame(data=values, index=rows).transpose()
    print(mut_table.to_csv(sep="\t"))

    return t_z


def get_rbd_statistics(table: pd.DataFrame):

    cols = [
        "SD1:mCSM-PPI2", "SD2:mCSM-PPI2",
        "SD1:MutaBind2", "SD2:MutaBind2",
        "SD1:SAAMBE-3D", "SD2:SAAMBE-3D",
    ]
    rows = list()
    values = list()

    wt_filtered = table.filter(like="_SWT:", axis=0)

    for swt_index in wt_filtered.index:

        # Create new row name.

        row_name = swt_index.replace("_SWT:", ":")

        # Calculate dddG values:
        # dddg_p[i]_sd[j] is the change in ddG predicted by the [i]th predictor (i = 1, 2, 3)
        #  when we move from the WT S-RBD to the Delta[j] variant S-RBD (j = 1, 2).

        sd1_index = swt_index.replace("_SWT:", "_SD1:")
        sd2_index = swt_index.replace("_SWT:", "_SD2:")

        if sd1_index not in table.index or sd2_index not in table.index:
            continue

        dddg_p1_sd1 = table["mCSM-PPI2"][sd1_index] - wt_filtered["mCSM-PPI2"][swt_index]
        dddg_p1_sd2 = table["mCSM-PPI2"][sd2_index] - wt_filtered["mCSM-PPI2"][swt_index]

        dddg_p2_sd1 = table["MutaBind2"][sd1_index] - wt_filtered["MutaBind2"][swt_index]
        dddg_p2_sd2 = table["MutaBind2"][sd2_index] - wt_filtered["MutaBind2"][swt_index]

        dddg_p3_sd1 = table["SAAMBE-3D"][sd1_index] - wt_filtered["SAAMBE-3D"][swt_index]
        dddg_p3_sd2 = table["SAAMBE-3D"][sd2_index] - wt_filtered["SAAMBE-3D"][swt_index]

        rows.append(row_name)
        values.append([dddg_p1_sd1, dddg_p1_sd2, dddg_p2_sd1, dddg_p2_sd2, dddg_p3_sd1, dddg_p3_sd2])

    dddg_table = pd.DataFrame(data=values, index=rows, columns=cols)

    # Print statistics
    for col_name in cols:

        current_table = dddg_table[col_name].sort_values()

        print(f"\nvariant:predictor = {col_name}")

        print("- HEAD -\n" + current_table.head(n=60).to_string())
        print("- TAIL -\n" + current_table.tail(n=60).to_string())
