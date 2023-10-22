from pathlib import Path
from typing import Tuple, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.neighbors import KernelDensity
from sklearn.metrics import roc_curve, auc
from scipy.stats import spearmanr


PLOT_DIR = Path("plots")
plt.rcParams.update({"font.size": 16})


def series_to_density(
        ser: Union[np.ndarray, pd.Series],
        min_value=None,
        max_value=None,
        n_of_point=None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Creates density datapoints with KDE from a series of data.

    :param ser: A 1D Pandas Series instance or a numpy ndarray instance.
    :param min_value: The minimum value for the x range.
    :param max_value: The maximum value for the x range.
    :param n_of_point: The number of points in the x range.
    :return: The x and y values of the datapoints.
    """

    if type(ser) is pd.Series:
        ser: np.ndarray = ser.to_numpy()

    min_value = np.min(ser) if min_value is None else min_value
    max_value = np.max(ser) if max_value is None else max_value
    n_of_points = 1000 if n_of_point is None else n_of_point

    x_values = np.arange(min_value, max_value, (max_value - min_value) / n_of_points)

    kde = KernelDensity(bandwidth="scott")
    kde.fit(ser[:, np.newaxis])
    y_values = kde.score_samples(x_values[:, np.newaxis])

    return x_values, np.exp(y_values)


def plot_distributions(table: pd.DataFrame) -> None:
    """
    Creates plots that help t visualize the distribution of the predicted ddG values.

    :param table:
    """

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


def categoric_mutations(ser: pd.Series, category_resi: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Filters the mutations inside a series based on whether they are "in cat -> not in cat" (out of category)
    or "not in cat -> cat" (into category) mutations.
    "in cat -> in cat" and "not in cat -> not in cat" mutations are thrown away.

    :param ser: The Pandas Series instance
    :param category_resi: The OLC of the amino acids inside the specific category.
    :return: The filtered values and the "out of category mutation" property (mask) of those values.
    """

    values = list()
    out_of_cat = list()
    for ser_idx in ser.index:

        mut = ser_idx.split(":")[1]
        wt_in_cat = mut[0] in category_resi
        mut_in_cat = mut[-1] in category_resi

        if wt_in_cat and not mut_in_cat:
            values.append(ser[ser_idx])
            out_of_cat.append(True)
        elif not wt_in_cat and mut_in_cat:
            values.append(ser[ser_idx])
            out_of_cat.append(False)

    return np.array(values, dtype=float), np.array(out_of_cat, dtype=bool)


def categoric_densities(values: np.ndarray, mask: np.ndarray) -> Tuple[np.ndarray, ...]:
    """
    Does the same as series_to_density, but the minimum and maximum values are coming from
    the full value set and densities are calculated for "out of category mutations" as well as
    for "into-category mutations".

    :param values: The predicted ddG values.
    :param mask: A boolean vector with values of True for "out of category mutations" and False for
     "into-category mutations".
    :return: The x values and the corresponding "out of category mutation" densities and
     "into-category mutation" densities.
    """

    min_value = np.min(values)
    max_value = np.max(values)

    x_values, y1_values = series_to_density(values[mask], min_value, max_value)
    _, y2_values = series_to_density(values[~mask], min_value, max_value)

    return x_values, y1_values, y2_values


def plot_correlation_scatter(table: pd.DataFrame) -> None:
    """
    Creates cross-scatter plots between the different predictors to view the correlation between them.

    :param table: The Pandas DataFrame instance containing the predicted ddG values.
    """

    fig, ax = plt.subplots(1, 3)
    fig.suptitle("Scatter plots between the different predictors to assess correlation")
    fig.subplots_adjust(wspace=0.3)
    fig.set_size_inches(13, 5)

    for axis in ax:
        axis.set_aspect("equal")
        axis.grid(visible=True, linestyle="dashed")
        axis.set_axisbelow(True)

    spr1 = spearmanr(table["mCSM-PPI2"].to_numpy(), table["MutaBind2"].to_numpy()).statistic
    spr2 = spearmanr(table["mCSM-PPI2"].to_numpy(), table["SAAMBE-3D"].to_numpy()).statistic
    spr3 = spearmanr(table["MutaBind2"].to_numpy(), table["SAAMBE-3D"].to_numpy()).statistic

    ax[0].set_title(f"SpR = {spr1:.3f}")
    ax[1].set_title(f"SpR = {spr2:.3f}")
    ax[2].set_title(f"SpR = {spr3:.3f}")

    ax[0].scatter(table["mCSM-PPI2"], table["MutaBind2"], color="black", marker=".", alpha=0.3, edgecolors=None)
    ax[1].scatter(table["mCSM-PPI2"], table["SAAMBE-3D"], color="black", marker=".", alpha=0.3, edgecolors=None)
    ax[2].scatter(table["MutaBind2"], table["SAAMBE-3D"], color="black", marker=".", alpha=0.3, edgecolors=None)

    dim_str = "$\\Delta\\Delta G$ / $\\mathrm{kcal} \\cdot \\mathrm{mol}^{-1}$"

    ax[0].set_xlabel(f"mCSM-PPI2 {dim_str}")
    ax[0].set_ylabel(f"MutaBind2 {dim_str}")

    ax[1].set_xlabel(f"mCSM-PPI2 {dim_str}")
    ax[1].set_ylabel(f"SAAMBE-3D {dim_str}")

    ax[2].set_xlabel(f"MutaBind2 {dim_str}")
    ax[2].set_ylabel(f"SAAMBE-3D {dim_str}")

    fig_name = "correlation_scatter.png"

    fig.savefig(PLOT_DIR / fig_name, dpi=300)


def plot_categoric_mutations(table: pd.DataFrame, category_resi: str):

    t1, cat_mask = categoric_mutations(table["mCSM-PPI2"], category_resi)
    t2, _ = categoric_mutations(table["MutaBind2"], category_resi)
    t3, _ = categoric_mutations(table["SAAMBE-3D"], category_resi)

    fig, ax = plt.subplots(1, 3)
    fig.suptitle(f"Out of and into {category_resi} group separation of predictors")
    fig.subplots_adjust(wspace=0.3, top=0.8)
    fig.set_size_inches(14, 6)

    x1_values, y1_out_values, y1_in_values = categoric_densities(t1, cat_mask)
    ax[0].set_title("mCSM-PPI2")
    ax[0].plot(x1_values, y1_out_values, color="blue", label="out of group")
    ax[0].plot(x1_values, y1_in_values, color="red", label="into group")

    x2_values, y2_out_values, y2_in_values = categoric_densities(t2, cat_mask)
    ax[1].set_title("MutaBind2")
    ax[1].plot(x2_values, y2_out_values, color="blue", label="out of group")
    ax[1].plot(x2_values, y2_in_values, color="red", label="into group")

    x3_values, y3_out_values, y3_in_values = categoric_densities(t3, cat_mask)
    ax[2].set_title("SAAMBE-3D")
    ax[2].plot(x3_values, y3_out_values, color="blue", label="out of group")
    ax[2].plot(x3_values, y3_in_values, color="red", label="into group")

    for axis in ax:

        axis.set_xlabel("$\\Delta\\Delta G$ / $\\mathrm{kcal} \\cdot \\mathrm{mol}^{-1}$")
        axis.set_ylabel("Density")
        axis.legend(loc="upper center")
        y_min, y_max = axis.get_ylim()
        axis.set_ylim(y_min, 1.3 * y_max)

    fig_name = f"categoric_separation_density_{category_resi}.png"

    fig.savefig(PLOT_DIR / fig_name, dpi=300)

    # Get ROC-AUC for the separation power of the given category

    plt.close(fig)
    fig, ax = plt.subplots(1, 3)
    fig.suptitle(
        f"ROC-AUC description of the classification ability of $\\Delta\\Delta G$ values\n"
        f"with respect to into- and out of {category_resi} group mutations"
    )
    fig.subplots_adjust(wspace=0.4)
    fig.set_size_inches(13, 5)

    ax[0].set_title("mCSM-PPI2")
    ax[1].set_title("MutaBind2")
    ax[2].set_title("SAAMBE-3D")

    tick_posi = np.arange(0, 1.2, 0.2)

    for axis in ax:
        axis.set_xlabel("false out of group rate"),
        axis.set_ylabel("true out of group rate"),
        axis.set_xticks(tick_posi),
        axis.set_yticks(tick_posi),
        axis.set_xlim(-0.05, 1.05),
        axis.set_ylim(-0.05, 1.05),
        axis.set_aspect("equal"),
        axis.grid()

    fpr1, tpr1, _ = roc_curve(np.array(cat_mask, dtype=int), t1)
    fpr2, tpr2, _ = roc_curve(np.array(cat_mask, dtype=int), t2)
    fpr3, tpr3, _ = roc_curve(np.array(cat_mask, dtype=int), t3)

    auc1 = auc(fpr1, tpr1)
    auc2 = auc(fpr2, tpr2)
    auc3 = auc(fpr3, tpr3)

    ax[0].plot(fpr1, tpr1, color="black")
    ax[1].plot(fpr2, tpr2, color="black")
    ax[2].plot(fpr3, tpr3, color="black")

    bbox_style = {"boxstyle": "round", "facecolor": "white", "alpha": 0.5}
    alignment_style = {"verticalalignment": "bottom", "horizontalalignment": "right"}
    ax[0].text(0.95, 0.05, f"{auc1:.3f}", transform=ax[0].transAxes, bbox=bbox_style, **alignment_style)
    ax[1].text(0.95, 0.05, f"{auc2:.3f}", transform=ax[1].transAxes, bbox=bbox_style, **alignment_style)
    ax[2].text(0.95, 0.05, f"{auc3:.3f}", transform=ax[2].transAxes, bbox=bbox_style, **alignment_style)

    fig_name = f"roc_auc_{category_resi}.png"

    fig.savefig(PLOT_DIR / fig_name, dpi=300)
