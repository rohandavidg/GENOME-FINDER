#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
setting up logging
"""

import logging
import sys
import time
import datetime


log_filename = sys.argv[0]

def main():
    logger = configure_logger(log_filename)


def configure_logger(log_filename):
    """
    setting up logging
    """
    logger = logging.getLogger(log_filename)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(time.strftime(log_filename+"-%Y%m%d.log"))
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s'\t'%(name)s'\t'%(levelname)s'\t'%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


if __name__ == "__main__":
    main()
