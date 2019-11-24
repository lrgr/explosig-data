import sys
import logging
import pandas as pd

FORMAT = '%(levelname)-10s: %(message)s'
FORMAT_WITH_TIME = '%(asctime)s ' + FORMAT

def get_logger(console_verbosity=logging.INFO, outfile=None):
    if outfile is not None:
        logging.basicConfig(
            level=logging.DEBUG,
            format=FORMAT_WITH_TIME, 
            filename=outfile,
            filemode='w'
        )
        console = logging.StreamHandler()
        console.setLevel(console_verbosity)
        formatter = logging.Formatter(FORMAT)
        console.setFormatter(formatter)

        logging.getLogger('').addHandler(console)
    else:
        logging.basicConfig(
            level=logging.DEBUG,
            format=FORMAT_WITH_TIME
        )

    return logging.getLogger(__name__)

def get_df_drop_message(col, reason, df_0, df_1):
    num_rows = df_0.shape[0] - df_1.shape[0]
    return "Dropping %i rows because %s in %s column" % (num_rows, reason, col)