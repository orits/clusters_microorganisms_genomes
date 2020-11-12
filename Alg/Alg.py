from Alg import functions


# this function will perform the algorithm stating.
# didnt return a value, just run it.
def alg(data, n_threads, completeness_threshold, contamination_threshold, genome_group_size,
        ani_threshold, linkage_threshold):
    # state 1 - check-m.
    functions.check_m(data, n_threads)

    # state 2 - init.
    genomes_list = functions.init_genomes(data, completeness_threshold, contamination_threshold)

    # state 3 - genomes_list sorted by the quality value.
    functions.sorted_genomes_list(genomes_list)

    # state 4 - clustering.
    clusters_list = functions.cluster_the_genomes(genomes_list, genome_group_size, n_threads, ani_threshold)

    # state 5 - sifting the clustering.
    clusters_list = functions.calculate_and_sift_link_of_clusters(clusters_list, n_threads, ani_threshold,
                                                                  linkage_threshold)

    # state 6 - representative.
    representatives_list = functions.select_representative_of_cluster(clusters_list)

    # state 7 - gtdbtk.
    functions.gtdbtk(clusters_list, representatives_list)

    # state 8 - create genome_cluster list.
    functions.create_final_file(clusters_list)
