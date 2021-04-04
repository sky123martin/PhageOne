import os
import multiprocessing 


class Config(object):
    # needed for CRSV which is used for forms
    SECRET_KEY = os.urandom(32) 

    HOMOLOGOUS_LENGTH = 100
    PROCESSES = multiprocessing.cpu_count()-2