"""
A class for defining an Cluster.
The c'tor receives the following arguments:
    genome: genome that create the cluster, temporary representative genome.
"""


class Cluster:
    # Static Class Attribute
    counter = 0

    # Initializer / Instance Attributes
    def __init__(self, genome):
        self.id = Cluster.counter
        self.genomes_list = []
        self.representative_genome = genome
        self.taxonomy = ""

        self.genomes_list.append(genome)
        Cluster.counter += 1

    # this function will select a representative to the cluster between all the genomes in the cluster.
    # all the genomes with the highest link score will be a candidates to be a representative.
    # from the candidates will chose one genome that  his median is the highest of all candidates.
    # return the genome that selected.
    def select_representative(self):
        cluster_genome_len = len(self.genomes_list)
        maximum_linked_val = -1
        if cluster_genome_len == 1:
            candidates_list = []
            for genome in self.genomes_list:
                if genome.linked > maximum_linked_val:
                    maximum_linked_val = genome.linked
                    candidates_list.clear()
                    candidates_list.append(genome)
                elif genome.linked == maximum_linked_val:
                    candidates_list.append(genome)

            maximum_ani_score = (-1, -1)
            for index, candidate in enumerate(candidates_list):
                if len(candidate.linked_set):
                    median_ani_score = candidate.clac_avg_ani_score()
                    if median_ani_score > maximum_ani_score[0]:
                        maximum_ani_score = (median_ani_score, index)

            self.representative_genome = candidates_list[maximum_ani_score[1]]
