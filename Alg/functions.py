import os
import sys
import logging
import itertools
from Alg.Genome import Genome
from Alg.Cluster import Cluster
import general_functions


# global var - logger.
logger = logging.getLogger("root")


# this function will check the quality of the genomes.
# didnt return a value, but will create an out file called "checkm.stdout".
def check_m(data, n_threads):
    logger.info("# state 1 #")
    try:
        if os.path.isdir("checkm"):
            os.rmdir("checkm")

        os.mkdir("checkm")
        os.mkdir("checkm/proteins")
        # create soft-link to all proteins files.

        for index, row in enumerate(data):
            is_valid, error_message = general_functions.check_is_legal_input_row(row)
            if not is_valid:
                logger.error("input row isn't valid (row index in file:" + str(index) + ") " + error_message)
                print("input file error, check log file.")
                sys.exit(0)

            path = row[2]
            name = row[0] + ".faa"

            os.system("ln -s " + path + " checkm/proteins/" + name)

    except OSError as e:
        line = sys.exc_info()[-1].tb_lineno
        logger.error(e.strerror + " (line: " + str(line) + ")")

    logger.info("# state 1")
    logger.info("starting checkm")
    os.system("~itaish/software/bin/checkm.lineage_wf.sh checkm/proteins checkm --ncpus=" + str(n_threads))
    logger.info("ending checkm")


# this function will initialize the genomes form the given data arg.
# didnt return a value.
def init_genomes(data, completeness_threshold, contamination_threshold):
    logger.info("# state 2 #")
    if os.path.isfile('checkm/checkm.stdout'):
        os.system("cat checkm/checkm.stdout  > checkm/checkm.stdout.txt")

    else:
        logger.error("file not exist")
        print("can't open file, exit from the program.")
        sys.exit(0)

    checkm_stdout_file = general_functions.read_file("checkm.stdout.txt")[3:-1]
    start_len = len(checkm_stdout_file)

    if start_len > 0:
        genomes_list = []
        logger.info("# state 2")
        logger.info("starting the mapping genomes process")
        for index, line in enumerate(checkm_stdout_file):
            completeness = float(line[12])
            contamination = float(line[13])
            if completeness >= float(completeness_threshold) and contamination <= float(contamination_threshold):
                name = data[index][0]
                path = data[index][1]
                genomes_list.append(Genome(name, path, completeness, contamination))

        removed_genomes = start_len - len(genomes_list)
        logger.info("ending the mapping genomes process - " + str(removed_genomes) + " genomes was removed")
        return genomes_list
    else:
        logger.error("the file: checkm.stdout.txt is empty!!")
        sys.exit(0)


# this function will sorted the genomes list by the genome quality
# in decreasing order.
# didnt return a value, just sort.
def sorted_genomes_list(genomes_list):
    logger.info("# state 3 #")
    sorted(genomes_list, key=lambda t: t.quality, reverse=True)


# this function will clustering all the given genomes_list to clusters.
# will cluster by taking a subgroup (from top) and run fastANI against the the bigger group the include the subgroup.
# return list of all created clusters.
def cluster_the_genomes(genomes_list, genome_group_size, n_threads, ani_threshold):
    logger.info("# state 4 #")
    clusters_list = []
    genome_cluster_map = list(map(lambda gen: (gen.path, False), genomes_list))
    genomes_to_remove = []
    logger.info("starting to clustering the genomes")

    try:
        if not os.path.isdir("init_clusters"):
            os.mkdir("init_clusters")

    except OSError as e:
        line = sys.exc_info()[-1].tb_lineno
        logger.error(e.strerror + " (line: " + str(line) + ")")

    need_clustering = True
    iteration = 0

    logger.info("starting to clustering the genomes")

    while need_clustering:
        path = "init_clusters/round_" + str(iteration)
        iteration += 1

        try:
            if not os.path.isdir(path):
                os.mkdir(path)

        except OSError as e:
            line = sys.exc_info()[-1].tb_lineno
            logger.error(e.strerror + " (line: " + str(line) + ")")

        if len(genomes_list) < genome_group_size:
            genome_group_size = len(genomes_list)

        # state a.

        # take the top genome_group_size of all genomes_list.
        queries_list = list(itertools.islice(genomes_list, genome_group_size))
        create_path_file_from_genomes(path + "/queries_list.txt", queries_list, False)

        # state b.
        create_path_file_from_genomes(path + "/subjects_list", genomes_list, False)

        # state c.
        outfile = path + "/fastAniOutFile.txt"
        os.system(
            "fastANI --ql " + path + "/queries.list.txt" + " --rl " + path + "/subjects.list.txt" + " -o " +
            outfile + " -t " + str(n_threads) + " > " + outfile + ".stderr 2> " + outfile + ".stdout")

        # state d.
        outfile = general_functions.read_file(outfile)

        if not outfile[-1]:  # remove last empty line from the output file.
            outfile = outfile[:-1]

        if outfile:  # file not empty.
            running_outfile = True
            while running_outfile:
                left_genome = outfile[0][0]
                is_left_cluster, left_index = check_if_genome_clustered(left_genome, genome_cluster_map)
                if not is_left_cluster:  # if left genome isn't clustered.
                    genome = genomes_list[left_index]  # get the left genome object by is index.
                    current_cluster = Cluster(genome)
                    genome_cluster_map[left_index] = (
                        genome_cluster_map[left_index][0], True)  # set the left genome as clustered.
                    if left_index >= genome_group_size:
                        genomes_to_remove.append(left_index)

                    if left_genome == outfile[0][1]:
                        outfile = outfile[1:]

                    sub_file, outfile = partition(lambda sub_row: sub_row[0] == left_genome and not (
                            sub_row[0] == left_genome and sub_row[1] == left_genome), outfile)

                    if not outfile:
                        running_outfile = False

                    if sub_file:
                        for index, sub_file_row in enumerate(sub_file):
                            right_genome = sub_file_row[1]
                            is_right_cluster, right_index = check_if_genome_clustered(right_genome, genome_cluster_map)
                            # if right genome isn't clustered, it may cluster to the current cluster.
                            if not is_right_cluster:
                                ani_score = float(sub_file_row[2])
                                if ani_score >= ani_threshold:
                                    current_cluster.genomes_list.append(genomes_list[right_index])
                                    genome_cluster_map[right_index] = (
                                        genome_cluster_map[right_index][0], True)  # set the right genome as clustered.
                                    if right_index >= genome_group_size:
                                        genomes_to_remove.append(right_index)

                    else:
                        logger.error("the sub file is empty!!")

                    clusters_list.append(current_cluster)  # add the cluster to the clusters list.
        else:
            logger.error("the file: " + path + "/fastAniOutFile.txt is empty!!")
            sys.exit(0)

        # state e.
        logger.info("# state 4 - e")
        genomes_list = update_genome_list(genomes_list, genomes_to_remove, genome_group_size)
        genome_cluster_map = update_genome_cluster_map(genome_cluster_map, genomes_to_remove, genome_group_size)
        genomes_to_remove = []  # set to new list for the next run.
        need_clustering = len(genomes_list) > 0

    return clusters_list


# this function will calculate the link score of all genomes and mapping the those that didnt compliance to conditions.
# will create a new folder to each cluster, then will run fastANI all X all,
# func: calc_links_in_cluster - will calculate the score,
# func: sifting_low_link_genome_in_cluster - will mapping all of those that didnt compliance to conditions.
def calculate_and_sift_link_of_clusters(clusters_list, n_threads, ani_threshold, linkage_threshold):
    singleton_cluster_list = []

    logger.info("# state 5 #")
    logger.info("starting to mapping the clusters")

    try:
        if not os.path.isdir("clusters"):
            os.mkdir("clusters")
    except OSError as e:
        line = sys.exc_info()[-1].tb_lineno
        logger.error(e.strerror + " (line: " + str(line) + ")")

    for cluster_index, cluster in enumerate(clusters_list):
        cluster_genome_len = len(cluster.genomes_list)
        path = "clusters/cluster_" + str(cluster.id)

        try:
            if not os.path.isdir(path):
                os.mkdir(path)

        except OSError as e:
            line = sys.exc_info()[-1].tb_lineno
            logger.error(e.strerror + " (line: " + str(line) + ")")

        genomes_path_list = path + "/genome_list_.before_maping.txt"
        create_path_file_from_genomes(genomes_path_list, cluster.genomes_list, False)

        if cluster_genome_len > 1:
            outfile = path + "/fastAniOutFile.txt"
            os.system(
                "fastANI --ql " + genomes_path_list + " --rl " + genomes_path_list + " -o " + outfile + " -t "
                + str(n_threads) + " --matrix" + " > " + outfile + ".stderr 2> " + outfile + ".stdout")

            outfile = general_functions.read_file(outfile + ".matrix")
            if not outfile[-1]:
                outfile = outfile[:-1]

            # calculate the linked score of all the genomes in the cluster.
            calc_links_in_cluster(cluster, outfile, ani_threshold, path)

            # create a friendship file before the mapping state.
            friendship_outfile(path + "/friendshipOutFileBeforeSift.txt", cluster.genomes_list)

            # sift the linked score of all the genomes in the cluster less the linkage_threshold.
            has_changes = sifting_low_link_genome_in_cluster(cluster, singleton_cluster_list, linkage_threshold)
            create_path_file_from_genomes(path + "/genomeListAfterSift.txt", cluster.genomes_list, False)

            if has_changes:  # create new and update cluster genomes list.
                friendship_outfile(path + "/friendshipOutFileAfterSift.txt", cluster.genomes_list)

    logger.info("ending to mapping the clusters - " + str(len(singleton_cluster_list))
                + " new singleton cluster was create")
    return singleton_cluster_list


# this function will select a representative to the cluster, by called method of class "Cluster"
# select_representative(),
# return a list of all representatives that selected.
def select_representative_of_cluster(clusters_list):
    logger.info("# state 6 #")
    logger.info("starting selecting clusters representatives")
    representatives_list = []
    for cluster in clusters_list:
        if len(cluster.genomes_list) > 1:
            cluster.select_representative()
        representatives_list.append(cluster.representative_genome)

    logger.info("ending selecting clusters representatives")
    return representatives_list


# this function will preform a gtdbtk state, that helper to classification the
# belonging of the representative genome to the tree of life
# didnt return a value, just set the cluster the taxonomy of the representative genome.
def gtdbtk(clusters_list, representatives_list):
    logger.info("# state 7 #")
    logger.info("create representatives_list file")

    if not os.path.isdir("gtdbtk"):
        os.mkdir("gtdbtk")

    create_path_file_from_genomes("gtdbtk/representatives_list.txt", representatives_list, True)

    out_dir = "gtdbtk"
    logger.info("running the gtdbtk comment")

    os.system("gtdbtk classify_wf --batchfile gtdbtk/representatives_list.txt --out_dir gtdbtk --cpus 20" + " > " +
              "gtdbtk/gtdbtk.stderr 2> " + "gtdbtk/gtdbtk.stdout")

    if os.path.isfile(out_dir + "/gtdbtk.ar122.summary.tsv"):
        outfile_one = general_functions.read_file(out_dir + "/gtdbtk.ar122.summary.tsv")[1:]
        logger.info("starting parse the the gtdbtk.ar122.summary.tsv file and create out file one.")
        for row_file_one in outfile_one:
            index = int(row_file_one[0])
            clusters_list[index].taxonomy = row_file_one[1]

        logger.info("ending parse the the gtdbtk.ar122.summary.tsv file")

    if os.path.isfile(out_dir + "/gtdbtk.bac120.summary.tsv"):
        outfile_two = general_functions.read_file(out_dir + "/gtdbtk.bac120.summary.tsv")[1:]
        logger.info("starting parse the the gtdbtk.bac120.summary.tsv file")
        for row_file_two in outfile_two:
            index = int(row_file_two[0])
            clusters_list[index].taxonomy = row_file_two[1]

        logger.info("ending parse the the gtdbtk.bac120.summary.tsv file")

    with open("representatives_taxonomy_list", 'w') as output:
        for cluster in clusters_list:
            output.write(str(cluster.id) + "\t" + cluster.representative_genome.name + "\t" + cluster.taxonomy + '\n')


# this function will create the first script output file.
# create a file with <genome-name>	<cluster-ID> for all genomes at all clusters.
# didnt return a value, just create a file.
def create_final_file(clusters_list):
    logger.info("# state 8 #")
    with open("genome_cluster_list", 'w') as output:
        for cluster in clusters_list:
            for genome in cluster.genomes_list:
                output.write(genome.name + "\t" + str(cluster.id) + "\n")

    logger.info("finish...")


# this function will create new file that contain only paths form genomes list.
# didnt return a value, create new file named "file_name".
def create_path_file_from_genomes(file_name, genomes_list, with_index):
    genomes_list_len = len(genomes_list)
    if genomes_list_len == 0:
        logger.error("genomes_list is empty!!")

    with open(file_name, 'w') as output:
        for index, genome in enumerate(genomes_list):
            if with_index:
                output.write(genome.path + '\t' + str(index))
            else:
                output.write(genome.path)

            if index < genomes_list_len - 1:
                output.write("\n")


# this function is a helper to the clustering state,
# search for tuple with the genome path and False value,
# return new tuple True/False and the index location,
# if didnt found return True and -1 index, else False, and founded index.
def check_if_genome_clustered(genome_path, genome_cluster_map):
    try:
        index = genome_cluster_map.index((genome_path, False))
        return False, index
    except ValueError:
        return True, -1


# this function will do partition at the given iterable, by the condition at the predicate.
# those element that compliance to the conditions stored at the trues list, else at the falses list.
# return trues and falses lists.
def partition(pred, iterable):
    trues = []
    falses = []
    for item in iterable:
        if pred(item):
            trues.append(item)
        else:
            falses.append(item)
    return trues, falses


# this function will update the genomes_list by get a list of genomes to remove.
# return a updated list of genomes without the genomes that all ready clustered.
def update_genome_list(genomes_list, genomes_to_remove, genome_group_size):
    new_genomes_list = genomes_list[genome_group_size:]
    for index in genomes_to_remove:
        del new_genomes_list[index - genome_group_size]

    return new_genomes_list


# this function will update the genome_cluster_map by get a list of genomes to remove.
# return a updated genome_cluster_map without the tuples that all ready used.
def update_genome_cluster_map(genome_cluster_map, genomes_to_remove, genome_group_size):
    new_genome_cluster_map = genome_cluster_map[genome_group_size:]
    for index in genomes_to_remove:
        del new_genome_cluster_map[index - genome_group_size]

    return new_genome_cluster_map


# this function will if indexs that send at the pred are friends,
# if yes, return tuple with the True , and the inner index at the genome linked_score list.
# if no, return a tuple with False, -1.
def check_are_we_friends(pred, iterable):
    for index, elem_tuple in enumerate(iterable):
        if pred(elem_tuple):
            return True, index
    return False, -1


# this function will remove the excess pairs of tuple from others that marked to delete.
# didnt return a value just update all friends of the index that removed that he is removed.
def delete_remaining_friends(index, cluster, index_set_to_remove):
    for inner_index, friend_tuple in enumerate(cluster.genomes_list[index].linked_score):
        friend_index = friend_tuple[0]
        if friend_index not in index_set_to_remove:
            are_we_friends, inner_index_friend = check_are_we_friends(
                lambda f_tuple: f_tuple[0] == index,
                cluster.genomes_list[friend_index].linked_score)
            if are_we_friends:
                cluster.genomes_list[friend_index].linked -= 1
                del cluster.genomes_list[friend_index].linked_score[inner_index_friend]


# this function will create a friendship outfile file that explanation all the friendships inside the cluster genomes.
# didnt return a value, just create a file.
def friendship_outfile(file_name, genomes_list):
    cluster_genome_len = len(genomes_list)
    with open(file_name, 'w') as output:
        for index, genome in enumerate(genomes_list):
            linked_percent_of_cluster = (genome.linked / (cluster_genome_len - 1)) * 100
            output.write(genome.path + "\t" + str(genome.linked) + "\t" + str(linked_percent_of_cluster) + "%\n")


# this function will calculate the score of genome in the cluster X all genomes at the cluster the if them will called
# friends, clac friendships by a given ani_score by the file outfile, (matrix lower bottom)
# if the score between genome A to genome B is "NA" so print an error to the log file the indicate a missing couple.
# didnt return a value, just update genomes: linked, linked_score, linked_set.
def calc_links_in_cluster(cluster, outfile, ani_threshold, path):
    grow_length_list = []
    with open(path + "/messing_ani.txt", 'w') as messing_outfile:
        for row_index, row in enumerate(outfile):
            if row_index == 0:
                continue

            elif row_index == 1:
                grow_length_list.append(cluster.genomes_list[row_index - 1])

            else:
                temp_width_list = row[1:]

                for col_index, col_val in enumerate(temp_width_list):
                    if col_val != "NA":
                        avg_ani_score = float(col_val)
                        if avg_ani_score >= ani_threshold:
                            cluster.genomes_list[row_index - 1].linked += 1
                            cluster.genomes_list[row_index - 1].linked_score.append((col_index, avg_ani_score))
                            cluster.genomes_list[row_index - 1].linked_set.add(col_index)

                            cluster.genomes_list[col_index].linked += 1
                            cluster.genomes_list[col_index].linked_score.append((row_index - 1, avg_ani_score))
                            cluster.genomes_list[col_index].linked_set.add(row_index - 1)
                    else:
                        # error log that couple are missing - to log file.
                        missing_left_genome = cluster.genomes_list[row_index - 1].path
                        missing_right_genome = cluster.genomes_list[col_index].path
                        messing_outfile.write("messing couple: " + str(missing_left_genome) + ",\t" +
                                              str(missing_right_genome))

    if not general_functions.is_non_zero_file(path + "/messing_ani.txt"):  # delete file if empty.
        os.remove(path + "/messing_ani.txt")


# this function will mapping the genomes inside the cluster, used a helper recursive method for collecting all the index
# that need to be delete because a index marked as weak link.
# return a True if change has make, or False is no genomes has mapped.
def sifting_low_link_genome_in_cluster(cluster, singleton_cluster_list, linkage_threshold):
    cluster_genome_len = len(cluster.genomes_list)
    index_set_to_remove = set()
    need_to_map = True

    while need_to_map:
        need_to_map = False
        for index, genome in enumerate(cluster.genome_list):  # search for week link percent.

            if index not in index_set_to_remove:  # avoid from removed index to re-checked again.
                temp_cluster_genome_len = cluster_genome_len - len(index_set_to_remove)

                if temp_cluster_genome_len > 1 and ((genome.linked / (
                        temp_cluster_genome_len - 1)) * 100) < linkage_threshold:  # genome need to sift.
                    index_set_to_remove.add(index)
                    temp_cluster_genome_len = cluster_genome_len - len(index_set_to_remove)

                    for index_value in genome.linked_set:
                        if temp_cluster_genome_len > 1 and ((cluster.genomes_list[index_value].linked / (
                                temp_cluster_genome_len - 1)) * 100) < linkage_threshold:  # genome need to sift.
                            need_to_map = True
                            break

    if len(index_set_to_remove) + 1 == cluster_genome_len:
        index_list_to_remove = list(range(1, cluster_genome_len))
        index_list_to_remove.reverse()

    else:
        index_list_to_remove = sorted(index_set_to_remove, reverse=True)

    if len(index_list_to_remove):
        for i, index_val in enumerate(index_list_to_remove):
            delete_remaining_friends(index_val, cluster, index_list_to_remove)
            cluster.genomes_list[index_val].linked = 0
            cluster.genomes_list[index_val].linked_score.clear()
            cluster.genomes_list[index_val].linked_set.clear()
            singleton_cluster_list.append(Cluster(cluster.genomes_list[index_val]))
            del cluster.genomes_list[index_val]

        if index_set_to_remove:
            return True
    return False
