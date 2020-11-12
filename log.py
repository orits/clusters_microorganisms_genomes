import logging


# this function create a singleton logger to help log the script.
# return the logger by name "root".
def setup_custom_logger(name):
    logging.basicConfig(filename="logfile.log", level=logging.DEBUG,
                        format="%(levelname)s %(asctime)s %(module)s - %(message)s",
                        datefmt="%d/%m/%Y %H:%M:%S", filemode='w')
    handler = logging.StreamHandler()
    logger = logging.getLogger(name)
    logger.addHandler(handler)
    return logger
