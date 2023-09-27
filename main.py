import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

from sklearn.neighbors import KernelDensity
from data_readers import read_mutabind2, read_saambe_3d, read_mcsm_ppi2

PLOT_DIR = Path("plots")


def get_statistics(table: pd.DataFrame):

    t_med = table.median()
    t_mad = (table - t_med).abs().median()
    t_z = (table - t_med) / t_mad
    t_z = t_z.mean(axis=1)
    t_z = t_z.sort_values()

    table_stats = pd.concat({
        "Minimum": (t_min := table.min()),
        "Maximum": (t_max := table.max()),
        "Range": t_max - t_min,
        "Mean": table.mean(),
        "StD": table.std(),
        "Median": (t_med := table.median()),
        "MAD": (table - t_med).abs().median()
    }, axis=1)

    print(table_stats, "\n")

    t_ranks = table.rank(axis=0)
    print(t_ranks)

    print("Pearson\'s correlation matrix:")
    print(table.corr(method="pearson"), "\n")

    print("Spearman\'s correlation matrix:")
    print(table.corr(method="spearman"), "\n")


def series_to_density(ser: pd.Series):
    """
    Creates density datapoints with KDE from a series of data.

    :param ser: The Pandas Series instance
    :return: The x and y values of the datapoints.
    """

    min_value = ser.min()
    max_value = ser.max()
    n_of_points = 1000
    x_values = np.arange(min_value, max_value, (max_value - min_value) / n_of_points)

    kde = KernelDensity(bandwidth="scott")
    kde.fit(ser.to_numpy()[:, np.newaxis])
    y_values = kde.score_samples(x_values[:, np.newaxis])

    return x_values, np.exp(y_values)


def categoric_mutations(ser: pd.Series, categories: str):
    """
    Filters the mutations inside a series based on whether they are "in cat -> not in cat" (green)
    or "not in cat -> cat" (blue) mutations. "in cat -> in cat" and "not in cat -> not in cat" mutations
    are thrown away. If the categories string is empty, it will simply return the values with the black color
    vector.

    :param ser: The Pandas Series instance
    :param categories: The OLC of the amino acids inside the specific category.
    :return: The filtered values and the colors associated with the mutation direction.
    """

    if categories == "":
        return ser.to_numpy(), np.array(["black" for _ in ser])

    colors = list()
    for mut in ser.index:

        mut = mut.split(":")[1]
        wt_in_cat = mut[0] in categories
        mut_in_cat = mut[-1] in categories

        if wt_in_cat and not mut_in_cat:
            colors.append("green")
        elif not wt_in_cat and mut_in_cat:
            colors.append("blue")
        else:
            colors.append(0)

    colors = np.array(colors, dtype=object)
    mask = colors != 0
    values = ser.to_numpy()[mask]
    colors = colors[mask]

    return values, colors


def plot_categoric_mutations(table: pd.DataFrame, categories: str):

    fig, ax = plt.subplots(1, 3)
    fig.suptitle(
        "Predicted $\\Delta\\Delta G$ values in $\\mathrm{kcal} \\cdot \\mathrm{mol}^{-1}$ dimension\n"
        "for the different predictors"
    )
    fig.subplots_adjust(wspace=0.3)
    fig.set_size_inches(13, 5)

    for axis in ax:
        axis.set_aspect("equal")
        axis.grid(visible=True, linestyle="dashed")
        axis.set_axisbelow(True)

    t1, colors = categoric_mutations(table["mCSM-PPI2"], categories)
    t2, _ = categoric_mutations(table["MutaBind2"], categories)
    t3, _ = categoric_mutations(table["SAAMBE-3D"], categories)

    ax[0].scatter(t1, t2, marker=".", c=colors, alpha=0.5, edgecolors=None)
    ax[1].scatter(t1, t3, marker=".", c=colors, alpha=0.5, edgecolors=None)
    ax[2].scatter(t2, t3, marker=".", c=colors, alpha=0.5, edgecolors=None)

    ax[0].set_xlabel("mCSM-PPI2")
    ax[0].set_ylabel("MutaBind2")

    ax[1].set_xlabel("mCSM-PPI2")
    ax[1].set_ylabel("SAAMBE-3D")

    ax[2].set_xlabel("MutaBind2")
    ax[2].set_ylabel("SAAMBE-3D")

    fig_name = "cross_scatter.png" if categories == "" else f"cross_scatter_{categories}.png"

    fig.savefig(PLOT_DIR / fig_name, dpi=300)


def plot_data(table: pd.DataFrame):

    plt.rcParams.update({"font.size": 14})

    # Create cross-scatter plots

    plot_categoric_mutations(table, "")
    plot_categoric_mutations(table, "FWYH")
    plot_categoric_mutations(table, "RHKDESTNQ")
    plot_categoric_mutations(table, "RK")
    plot_categoric_mutations(table, "DE")
    plot_categoric_mutations(table, "P")
    plot_categoric_mutations(table, "G")

    # Create boxplot

    fig, ax = plt.subplots()

    fig.suptitle("Boxplot for the different predictors")
    fig.set_size_inches(7, 7)

    # Help links:
    # https://matplotlib.org/3.1.1/gallery/statistics/boxplot.html
    # https://matplotlib.org/stable/gallery/statistics/boxplot_color.html#sphx-glr-gallery-statistics-boxplot-color-py
    # https://matplotlib.org/stable/gallery/statistics/boxplot_demo.html#sphx-glr-gallery-statistics-boxplot-demo-py
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.boxplot.html

    ax.boxplot(
        [table["mCSM-PPI2"], table["MutaBind2"], table["SAAMBE-3D"]],
        patch_artist=True,
        labels=["mCSM-PPI2", "MutaBind2", "SAAMBE-3D"],
        flierprops={
            "marker": ".", "markerfacecolor": "black", "markeredgewidth": 0, "alpha": 0.2
        },
        boxprops={
            "linewidth": 2.5, "facecolor": "cyan"
        },
        medianprops={
            "linewidth": 2.5, "color": "black"
        }
    )
    ax.set_ylabel("predicted $\\Delta\\Delta G$ / $\\mathrm{kcal} \\cdot \\mathrm{mol}^{-1}$")

    fig.savefig(PLOT_DIR / "boxplots.png", dpi=300)

    # Create density plots

    plt.close(fig)
    fig, ax = plt.subplots()

    fig.suptitle("Density plots for the different predictors")
    fig.set_size_inches(7, 7)

    x1_values, y1_values = series_to_density(table["mCSM-PPI2"])
    x2_values, y2_values = series_to_density(table["MutaBind2"])
    x3_values, y3_values = series_to_density(table["SAAMBE-3D"])

    ax.plot(x1_values, y1_values, label="mCSM-PPI2", c="red")
    ax.plot(x2_values, y2_values, label="MutaBind2", c="green")
    ax.plot(x3_values, y3_values, label="SAAMBE-3D", c="blue")

    ax.legend(loc="upper right")
    ax.set_xlabel("predicted $\\Delta\\Delta G$ / $\\mathrm{kcal} \\cdot \\mathrm{mol}^{-1}$")
    ax.set_ylabel("density")

    fig.savefig(PLOT_DIR / "densities.png", dpi=300)


def main():

    table = read_mcsm_ppi2()
    table = table.join(read_saambe_3d(), how="inner")
    table = table.join(read_mutabind2(), how="inner")

    # get_statistics(table)

    plot_data(table)


if __name__ == "__main__":
    main()
