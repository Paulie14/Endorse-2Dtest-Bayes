import os
import sys
import ruamel.yaml as yaml
import numpy as np

import flow_wrapper
from measured_data import MeasuredData
from mesh_factory import MeshFactory


def preprocess(config_dict):
    # prepare measured data as observations
    md = MeasuredData(config_dict)
    md.initialize()

    md.plot_all_data()
    md.plot_interp_data()

    conf_bayes = config_dict["surrDAMH_parameters"]

    boreholes = conf_bayes["observe_points"]
    times, values = md.generate_measured_samples(boreholes)

    config_bayes_file = config_dict["bayes_config_file"]
    yaml_handler = yaml.YAML()
    with open(config_bayes_file) as f:
        file_content = f.read()
    conf = yaml_handler.load(file_content)
    # print(conf.ca)

    npar = len(conf_bayes["parameters"])
    conf["no_parameters"] = npar
    conf["problem_parameters"]["noise_std"] = [100] # * len(values)
    conf["problem_parameters"]["observations"] = np.array(values).tolist()
    conf["problem_parameters"]["prior_mean"] = [0.0] * npar
    conf["problem_parameters"]["prior_std"] = [1.0] * npar
    conf["no_observations"] = len(values)
    conf["noise_type"] = "Gaussian_process"
    conf["noise_grid"] = np.array(times).tolist()
    conf["noise_parameters"] = [[30, 100]] * len(boreholes)
    conf["solver_module_path"] = os.path.join(config_dict["script_dir"], "flow_wrapper.py")
    conf["transformations"] = conf_bayes["parameters"]

    for i, par in enumerate(conf_bayes["parameters"]):
        if par["type"] is None:
            conf["problem_parameters"]["prior_mean"][i] = par["options"]["mu"]
            conf["problem_parameters"]["prior_std"][i] = par["options"]["sigma"]

    conf["no_solvers"] = int(np.round(0.5*(config_dict["metacentrum"]["chunks"] * config_dict["metacentrum"]["ncpus_per_chunk"]-2)))

    with open(config_bayes_file, 'w') as f:
        yaml_handler.dump(conf, f)

    MeshFactory.prepare_mesh(config_dict["geometry"], config_dict["common_files_dir"], cut_tunnel=True)


if __name__ == "__main__":

    output_dir = None
    len_argv = len(sys.argv)
    assert len_argv > 1, "Specify output directory!"
    if len_argv > 1:
        output_dir = sys.argv[1]

    config_dict = flow_wrapper.setup_config()
    preprocess(config_dict)


