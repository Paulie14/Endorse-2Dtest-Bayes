import os
import sys
import shutil
import ruamel.yaml as yaml

from flow123d_simulation import endorse_2Dtest


if __name__ == "__main__":
    with open(os.path.join(os.getcwd(), "config.yaml"), "r") as f:
        config_dict = yaml.safe_load(f)
    # self.config_dict["config_pbs"] = os.path.join(os.getcwd(), "config_PBS.yaml")

    work_dir = os.getcwd() + "/flow123d_sim"
    # Create working directory if necessary
    os.makedirs(work_dir, mode=0o775, exist_ok=True)

    config_dict["work_dir"] = work_dir
    config_dict["script_dir"] = os.getcwd()
    config_dict["_aux_flow_path"] = config_dict["local"]["flow_executable"].copy()
    config_dict["_aux_gmsh_path"] = config_dict["local"]["gmsh_executable"].copy()

    sim = endorse_2Dtest(config_dict, clean=True)
    sim.set_parameters(data_par=[6e-15, 0.17])
    res = sim.get_observations()

