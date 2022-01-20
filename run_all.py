import os
import sys
import shutil
import ruamel.yaml as yaml

import aux_functions


def setup():
    # create and cd workdir
    rep_dir = os.path.dirname(os.path.abspath(__file__))
    work_dir = os.path.join(rep_dir, "flow123d_sim")

    # Create working directory if necessary
    os.makedirs(work_dir, mode=0o775, exist_ok=True)
    os.chdir(work_dir)

    # read config file and setup paths
    with open(os.path.join(rep_dir, "config.yaml"), "r") as f:
        config_dict = yaml.safe_load(f)

    config_dict["work_dir"] = work_dir
    config_dict["script_dir"] = rep_dir

    clean = config_dict["clean_sample_dir"]

    # Files in the directory are used by each simulation at that level
    common_files_dir = os.path.join(work_dir, "common_files")
    aux_functions.force_mkdir(common_files_dir, force=clean)
    # copy common files
    for f in config_dict["copy_files"]:
        filepath = os.path.join(common_files_dir, f)
        if not os.path.isfile(filepath):
            shutil.copyfile(os.path.join(rep_dir, f), filepath)

    return config_dict


if __name__ == "__main__":

    # setup paths and directories
    config_dict = setup()

    # RUN THE MCMC SIMULATION
    # default parameters
    N = 2  # number of sampling processes
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

    # run sampling
    command = None
    if visualize:
        # os.error("Visualization not implemented.")
        # if os.path.exists("examples/visualization/" + problem_name + ".py"):
        #     command = "python3 examples/visualization/" + problem_name + ".py " + str(N)
        # else:
        #     command = "python3 examples/visualization/general_visualization.py " + str(N) + " " + problem_name
        command = "python3 examples/visualization/general_visualization.py " + str(N) + " " + problem_path
    else:
        if oversubscribe:
            opt = " --oversubscribe "
        else:
            opt = " "
        sampler = "python3 -m mpi4py surrDAMH/process_SAMPLER.py "
        solver = "python3 -m mpi4py surrDAMH/process_SOLVER.py " + problem_path + " "
        collector = "python3 -m mpi4py surrDAMH/process_COLLECTOR.py "

        if not config_dict["run_on_metacentrum"]:
            # possibly change to particular mpirun for testing
            # mpirun = "/usr/local/mpich_3.4.2/bin/mpirun"
            mpirun = "mpirun"
            command = mpirun + " -n " + str(N) + opt + sampler \
                      + ":" + " -n 1" + opt + solver + ":" + " -n 1" + opt + collector
        else:
            met = config_dict["metacentrum"]
            ncpus = met["ncpus"]
            lines = [
                '#!/bin/bash',
                '#PBS -S /bin/bash',
                '#PBS -l select=1:ncpus=' + str(ncpus+2) + ':mem=12gb',
                '#PBS -l walltime=' + met["walltime"],
                '#PBS -q ' + met["queue"],
                '#PBS -N ' + met["name"],
                '#PBS -j oe',
                '\n',
                'cd "' + config_dict["script_dir"] + '"',
                'image=$(./endorse_fterm image)',
                'sing_command="singularity exec -B ' + met["workspace_rel"] + '/:/' + met["workspace_rel"]
                        + ' docker://$image"',
                '\n',
                'bash_py="bash -c \'source ./venv/bin/activate &&"',
                'sampler="$bash_py ' + sampler + '\'"',
                'solver="$bash_py ' + solver + '\'"',
                'collector="$bash_py ' + collector + '\'"',
                '\n',
                'module load mpich-3.0.2-gcc',
                'which mpirun',
                'mpirun --version',
                'command="mpirun '
                        + '-n ' + str(ncpus) + ' $sing_command $sampler : '
                        + '-n 1 $sing_command $solver : '
                        + '-n 1 $sing_command $collector"',
                'echo $command',
                'eval $command'
            ]
            with open("shell_process.sh", 'w') as f:
                f.write('\n'.join(lines))

            lines = [
                '#!/bin/bash',
                'cd "' + config_dict["script_dir"] + '"',
                'image=$(./endorse_fterm image)',
                'sing_command="singularity exec -B ' + met["workspace_rel"] + '/:/' + met["workspace_rel"]
                + ' docker://$image"',
                '\n',
                'bash_py="bash -c \'source ./venv/bin/activate &&"',
                'sampler="$sing_command $bash_py ' + ' python3 preprocess.py ' + problem_path + '\'"',
                '\n',
                'echo $command',
                'eval $command'
            ]
            with open("shell_preprocess.sh", 'w') as f:
                f.write('\n'.join(lines))

    if config_dict["run_on_metacentrum"]:
        # TODO: run preprocess through singularity
        os.system(os.path.join(config_dict["work_dir"], "shell_preprocess.sh"))
        os.system("qsub " + os.path.join(config_dict["work_dir"], "shell_process.sh"))
    else:
        from preprocess import preprocess
        preprocess(config_dict, problem_path)

        rep_dir = config_dict["script_dir"]
        os.chdir(os.path.join(rep_dir, "surrDAMH"))
        print(command)
        # exit(0)
        os.system(command)
