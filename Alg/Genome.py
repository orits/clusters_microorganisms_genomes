import statistics

"""
A class for defining an Genome.
The c'tor receives the following arguments:
    name: genome identifier name from input file.
    path: genome fasta file path form input file.
    completeness: percent of completeness in the genome protein.
    contamination: percent of contamination in the genome protein.
"""


class Genome:

    # Initializer / Instance Attributes
    def __init__(self, name, path, completeness, contamination):
        self.name = name
        self.path = path
        self.completeness = completeness
        self.contamination = contamination
        self.quality = completeness - 5 * contamination
        self.linked = 0
        self.linked_score = []
        self.linked_set = set()

    # this function will clac the median ani score of a genome with all his friends.
    # return a median value all ani score with all his friends.
    def clac_median_ani_score(self):
        if len(self.linked_score) == 0:
            return 0

        temp, list_ani = zip(*self.linked_score)
        return statistics.median(list_ani)
