# -*- coding: utf-8 -*-



import os
import shutil
import ruamel.yaml as yaml

from flow123d_simulation import endorse_2Dtest


def force_mkdir(path, force=False):
    """
    Make directory 'path' with all parents,
    remove the leaf dir recursively if it already exists.
    :param path: path to directory
    :param force: if dir already exists then remove it and create new one
    :return: None
    """
    if force:
        if os.path.isdir(path):
            shutil.rmtree(path)
    os.makedirs(path, mode=0o775, exist_ok=True)


class Wrapper:
    def __init__(self, solver_id = 0):

        rep_dir = os.path.dirname(os.path.abspath(__file__))
        work_dir = os.path.join(rep_dir, "flow123d_sim")
        # Create working directory if necessary
        os.makedirs(work_dir, mode=0o775, exist_ok=True)
        os.chdir(work_dir)

        # read config file and setup paths
        with open(os.path.join(rep_dir, "config.yaml"), "r") as f:
            config_dict = yaml.safe_load(f)

        clean = config_dict["clean_sample_dir"]

        # Files in the directory are used by each simulation at that level
        common_files_dir = os.path.join(work_dir, "common_files")
        force_mkdir(common_files_dir, force=clean)
        config_dict["common_files_dir"] = common_files_dir

        # copy common files
        for f in config_dict["copy_files"]:
            filepath = os.path.join(common_files_dir, f)
            if not os.path.isfile(filepath):
                shutil.copyfile(os.path.join(rep_dir, f), filepath)

        config_dict["solver_id"] = solver_id
        config_dict["work_dir"] = work_dir
        config_dict["script_dir"] = rep_dir
        config_dict["_aux_flow_path"] = config_dict["local"]["flow_executable"].copy()
        config_dict["_aux_gmsh_path"] = config_dict["local"]["gmsh_executable"].copy()

        self.sim = endorse_2Dtest(config_dict, clean=clean)
        self.no_parameters = 2
        
    def set_parameters(self,data_par):
        self.sim.set_parameters(data_par)
        
    def get_observations(self):
        res = self.sim.get_observations()
        return res