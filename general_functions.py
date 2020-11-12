import os
import sys
import log
from datetime import datetime
from Alg import Alg

# create and configure the logger - starting to log...
logger = log.setup_custom_logger("root")


# this function will read a file by given path + filename.
# return the file data.
def read_file(file_name):
    if not os.path.isfile(file_name):
        logger.error("file: " + file_name + " not exist!!")
        sys.exit(0)

    data = []
    with open(file_name) as file:
        for line in file:
            data.append(line.split())

    return data


# this function run all the alg with his needed args.
# didnt return a value, just run.
def run(data, result_folder_path, n_threads, completeness_threshold, contamination_threshold, genome_group_size,
        ani_threshold, linkage_threshold):
    now = datetime.now()
    date_time = now.strftime("%d-%m-%Y_%H%M%S")
    result_folder_path = result_folder_path + "/result_" + str(date_time)
    try:
        if not os.path.isdir(result_folder_path):
            os.mkdir(result_folder_path)

        os.chdir(result_folder_path)

    except OSError as e:
        logger.error(e)
        print("please check your path: " + result_folder_path)
        sys.exit(1)

    logger.info("args:")
    logger.info("result_folder_path: " + str(result_folder_path))
    logger.info("n_threads: " + str(n_threads))
    logger.info("completeness_threshold: " + str(completeness_threshold))
    logger.info("contamination_threshold: " + str(contamination_threshold))
    logger.info("genome_group_size: " + str(genome_group_size))
    logger.info("ani_threshold: " + str(ani_threshold))
    logger.info("linkage_threshold: " + str(linkage_threshold))
    logger.info("starting the Alg run:")
    Alg.alg(data, n_threads, completeness_threshold, contamination_threshold, genome_group_size,
            ani_threshold, linkage_threshold)
    logger.info("ending the Alg run - take look at the result folder!")


# this function will check if path is a file and file isn't empty.
# return a boolean as answer.
def is_non_zero_file(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


# this function will check line from the input file of the program
# return boolean and a message.
def check_is_legal_input_row(input_row):
    n = len(input_row)
    if n != 3:
        return False, "legal row is 3 columns only!"

    if input_row[1][-4] == ".faa" or input_row[2][-4] == ".fna":
        return False, "legal row is tree columns like (genome_name tab path.fna tab path.faa)"

    return True, ""
