import os
import sys
import csv
import time

import flow_wrapper
from measured_data import MeasuredData

from run_all import setup
from preprocess import preprocess

from surrDAMH.surrDAMH.modules import visualization_and_analysis as va


def just_run_flow123d(measured_data, params_in, output_dir_in, boreholes_in):

    wrap = flow_wrapper.Wrapper(solver_id=0, output_dir=output_dir_in)

    for idx, pars in enumerate(params_in):
        wrap.set_parameters(data_par=pars)
        t = time.time()
        res, obs_data = wrap.get_observations()
        print("Flow123d res: ", res)
        if res >= 0:
            print(obs_data)
            measured_data.plot_comparison(obs_data, wrap.sim.sample_dir, boreholes_in)

        print("LEN:", len(obs_data))
        print("TIME:", time.time() - t)
        # if idx == 1:
        #     exit(0)


def get_best_accepted_params(config_dict_in, output_dir_in, count):

    data_dir = os.path.join(output_dir_in, "saved_samples", "config_mcmc_bayes")
    no_parameters = len(config_dict_in['surrDAMH_parameters']['parameters'])

    samples = va.Samples()
    samples.load_MH(data_dir, no_parameters)
    samples.calculate_properties()
    samples.load_MH_with_posterior(data_dir, no_parameters)

    # print(samples.x)
    # print(samples.x_compress)
    # print(samples.posteriors)

    # modus = samples.find_modus()
    # print("MODUS:\n", modus)

    modi = samples.find_n_modi(count)
    print(count, " best MODI:\n", modi)

    param_file = os.path.join(output_dir_in, "best_accepted_params.csv")
    params = []
    with open(param_file, 'w') as file:
        for post, i, li, mod in modi:
            params.append(mod)
            line = ','.join([str(s) for s in mod])
            file.write(line + "\n")

    return params


if __name__ == "__main__":

    # RUN THE MCMC SIMULATION
    # default parameters
    output_dir = None
    csv_data = None
    n_best_params = 0

    len_argv = len(sys.argv)
    assert len_argv > 2, "Specify output dir and parameters in csv file!"
    if len_argv > 1:
        output_dir = os.path.abspath(sys.argv[1])
    if len_argv > 2:
        file = sys.argv[2]
        if os.path.exists(file):
            csv_data = os.path.abspath(sys.argv[2])
        else:
            n_best_params = int(sys.argv[2])

    # setup paths and directories
    config_dict = setup(output_dir, can_overwrite=True)

    preprocess(config_dict)

    # prepare measured data as observations
    md = MeasuredData(config_dict)
    md.initialize()

    md.plot_all_data()
    md.plot_interp_data()

    boreholes = config_dict["surrDAMH_parameters"]["observe_points"]
    times, values = md.generate_measured_samples(boreholes)

    if csv_data:
        print("Reading parameters from CSV: ", csv_data)
        with open(csv_data, newline='') as csvfile:
            parameters = list(csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC))
    else:
        print("Getting " + str(n_best_params) + " best parameters.")
        parameters = get_best_accepted_params(config_dict, output_dir, n_best_params)

    print(parameters)

    print(boreholes)
    # JUST RUN FLOW123D FOR TESTING
    just_run_flow123d(md, parameters, output_dir, boreholes)
