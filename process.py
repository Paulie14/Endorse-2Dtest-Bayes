import os
import sys
import json
import time


def just_run_flow123d():
    from flow_wrapper import Wrapper
    wrap = Wrapper(solver_id=1)
    wrap.set_parameters(data_par=[6e-15, 0.17])
    t = time.time()
    res = wrap.get_observations()
    print(res)
    print("LEN:", len(res[1]))
    print("TIME:", time.time() - t)
    exit(0)


if __name__ == "__main__":

    # JUST RUN FLOW123D FOR TESTING
    # just_run_flow123d()

    # RUN THE MCMC SIMULATION
    # default parameters
    N = 4  # number of sampling processes
    oversubscribe = False  # if there are not enough slots available
    visualize = False  # True = only visualization
    problem_path = None

    len_argv = len(sys.argv)
    assert len_argv > 1, "Specify configuration json file!"
    if len_argv > 1:
        problem_path = sys.argv[1]
    if len_argv > 2:
        N = int(sys.argv[2])  # number of MH/DAMH chains
    if len_argv > 3:
        oversubscribe = sys.argv[3] == "oversubscribe"
        visualize = sys.argv[3] == "visualize"

    basename = os.path.basename(problem_path)
    problem_name, fext = os.path.splitext(basename)

    with open(problem_path) as f:
        conf = json.load(f)

    command = None
    if visualize:
        os.error("Visualization not implemented.")
        # if os.path.exists("examples/visualization/" + problem_name + ".py"):
        #     command = "python3 examples/visualization/" + problem_name + ".py " + str(N)
        # else:
        #     command = "python3 examples/visualization/general_visualization.py " + str(N) + " " + problem_name
    else:
        if oversubscribe:
            opt = " --oversubscribe "
        else:
            opt = " "
        sampler = " -n " + str(N) + opt + "python3 -m mpi4py surrDAMH/process_SAMPLER.py "
        solver = " -n 1" + opt + "python3 -m mpi4py surrDAMH/process_SOLVER.py " + problem_path + " "
        collector = " -n 1" + opt + "python3 -m mpi4py surrDAMH/process_COLLECTOR.py "
        if "surrogate_type" in conf.keys():
            command = "mpirun" + sampler + ":" + solver + ":" + collector
        else:
            command = "mpirun" + sampler + ":" + solver

    rep_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(os.path.join(rep_dir, "MCMC-Bayes-python"))
    print(command)
    os.system(command)
