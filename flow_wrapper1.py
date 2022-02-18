# -*- coding: utf-8 -*-


import os
import sys
import shutil
import ruamel.yaml as yaml
#import numpy as np

rep_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(rep_dir)

#import transformations as trans
from flow123d_simulation import endorse_2Dtest
import aux_functions


def setup_config():
    # create and cd workdir
    work_dir = os.path.join(rep_dir, "flow123d_sim")

    # read config file and setup paths
    with open(os.path.join(rep_dir, "config.yaml"), "r") as f:
        config_dict = yaml.safe_load(f)

    config_dict["work_dir"] = work_dir
    config_dict["script_dir"] = rep_dir
    config_dict["_aux_flow_path"] = config_dict["local"]["flow_executable"].copy()
    config_dict["_aux_gmsh_path"] = config_dict["local"]["gmsh_executable"].copy()

    # Files in the directory are used by each simulation at that level
    common_files_dir = os.path.join(work_dir, "common_files")
    config_dict["common_files_dir"] = common_files_dir

    return config_dict


def setup_dirs(config_dict):
    work_dir = config_dict["work_dir"]
    # Create working directory if necessary
    os.makedirs(work_dir, mode=0o775, exist_ok=True)
    os.chdir(work_dir)

    clean = config_dict["clean_sample_dir"]
    common_files_dir = config_dict["common_files_dir"]

    aux_functions.force_mkdir(common_files_dir, force=clean)
    # copy common files
    for f in config_dict["copy_files"]:
        filepath = os.path.join(common_files_dir, f)
        if not os.path.isfile(filepath):
            shutil.copyfile(os.path.join(rep_dir, f), filepath)


class Wrapper:
    def __init__(self, solver_id=0):

        config_dict = setup_config()

        config_dict["solver_id"] = solver_id

        clean = config_dict["clean_sample_dir"]
        self.sim = endorse_2Dtest(config_dict, clean=clean)
        
    def set_parameters(self, data_par):
        # conductivity = trans.normal_to_lognormal(data_par[0])
        # biot = trans.normal_to_beta(data_par[1], alfa=5, beta=5)
        # self.sim.set_parameters(np.array([conductivity, biot]))
        self.sim.set_parameters(data_par)
        
    def get_observations(self):
        res = self.sim.get_observations()
        return res
