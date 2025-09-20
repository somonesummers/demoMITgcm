# A cleaner version of run_utils.py

import os
import MITgcmutils # see https://mitgcm.readthedocs.io/en/latest/utilities/utilities.html
import re


def getLastIter(results_path):

    
#     if len(results_path)==0:
#         results_path = os.path.join(exp_dir, exp_name, 'results')
    
    # get most recent pickup file
    results_files = os.listdir(results_path)
    r = re.compile(r"^(pickup)\.\d+\.data")
    pickup_files = sorted(list(filter(r.search, results_files)))
    
    # get latest iter number
    r2 = re.compile(r"\d+")
    iterN_str = r2.search(pickup_files[-1]).group()
    
    return iterN_str