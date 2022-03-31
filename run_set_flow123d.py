import os
import sys
import csv
import time

import flow_wrapper
from measured_data import MeasuredData

from run_all import setup


def just_run_flow123d(measured_data, parameters, output_dir, boreholes):

    wrap = flow_wrapper.Wrapper(solver_id=0, output_dir=output_dir)

    for idx, pars in enumerate(parameters):
        wrap.set_parameters(data_par=pars)
        t = time.time()
        res, obs_data = wrap.get_observations()
        print("Flow123d res: ", res)
        if res >= 0:
            print(obs_data)
            measured_data.plot_comparison(obs_data, wrap.sim.sample_dir, boreholes)

        print("LEN:", len(obs_data))
        print("TIME:", time.time() - t)
        # if idx == 1:
        #     exit(0)


if __name__ == "__main__":

    # RUN THE MCMC SIMULATION
    # default parameters
    output_dir = None
    csv_data = None

    len_argv = len(sys.argv)
    assert len_argv > 2, "Specify output dir and parameters in csv file!"
    if len_argv > 1:
        output_dir = sys.argv[1]
    if len_argv > 2:
        csv_data = os.path.abspath(sys.argv[2])

    # setup paths and directories
    config_dict = setup(output_dir)

    # prepare measured data as observations
    md = MeasuredData(config_dict)
    md.initialize()

    md.plot_all_data()
    md.plot_interp_data()

    boreholes = config_dict["surrDAMH_parameters"]["observe_points"]
    times, values = md.generate_measured_samples(boreholes)

    with open(csv_data, newline='') as csvfile:
        parameters = list(csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC))

    print(boreholes)
    print(parameters)

    # JUST RUN FLOW123D FOR TESTING
    just_run_flow123d(md, parameters, output_dir, boreholes)
