import os
import sys
import json

import flow_wrapper
from measured_data import MeasuredData


def preprocess(config_dict, problem_path):
    # prepare measured data as observations
    md = MeasuredData(config_dict)
    md.initialize()

    md.plot_all_data()
    md.plot_interp_data()

    boreholes = config_dict["mcmc_observe_point"]
    times, values = md.generate_measured_samples(boreholes)

    with open(problem_path) as f:
        conf = json.load(f)
        conf["problem_parameters"]["observations"] = values
        conf["no_observations"] = len(values)
        conf["noise_type"] = "Gaussian_process"
        conf["noise_grid"] = times
        conf["noise_parameters"] = [[30, 50]] * 4
        conf["solver_module_path"] = os.path.join(config_dict["script_dir"], "flow_wrapper.py")

    with open(problem_path, "w") as f:
        json.dump(conf, f, indent=4)

    # TODO: move the mesh preparation here


if __name__ == "__main__":

    problem_path = None
    len_argv = len(sys.argv)
    assert len_argv > 1, "Specify configuration json file!"
    if len_argv > 1:
        problem_path = sys.argv[1]
    if len_argv > 2:
        output_dir = sys.argv[2]

    config_dict = flow_wrapper.setup_config()
    preprocess(config_dict, problem_path)


