from data_readers import read_mutabind2, read_saambe_3d, read_mcsm_ppi2
from plotters import plot_distributions, plot_correlation_scatter, plot_categoric_mutations
from statistics import get_global_statistics, get_rbd_statistics, get_z_statistics
from b_factor_coloring import modify_b_factor


def main():

    # Data loading

    table = read_mcsm_ppi2()
    table = table.join(read_saambe_3d(), how="inner")
    table = table.join(read_mutabind2(), how="inner")

    # Statistics calculation

    z_scores = get_z_statistics(table)
    get_rbd_statistics(table)
    get_global_statistics(table)

    # PDB b-factor modification based on z-scores

    modify_b_factor(z_scores)

    # Plotting

    plot_distributions(table)
    plot_correlation_scatter(table)

    plot_categoric_mutations(table, "FWYH")
    plot_categoric_mutations(table, "RHKDESTNQ")
    plot_categoric_mutations(table, "RK")
    plot_categoric_mutations(table, "DE")
    plot_categoric_mutations(table, "P")
    plot_categoric_mutations(table, "G")


if __name__ == "__main__":
    main()
