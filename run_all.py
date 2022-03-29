import os
import sys
import shutil
import ruamel.yaml as yaml

import aux_functions
from preprocess import preprocess
from surrDAMH.surrDAMH.configuration import Configuration


# this script is supposed to be dependent only on python packages present on any machine
# all other python scripts are later run inside docker container

def setup(output_dir):
    # create and cd workdir
    rep_dir = os.path.dirname(os.path.abspath(__file__))
    work_dir = os.path.abspath(output_dir)

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
        # if not os.path.isfile(filepath):
        shutil.copyfile(os.path.join(rep_dir, f), filepath)

    config_dict["bayes_config_file"] = os.path.join(common_files_dir,
                                                    config_dict["surrDAMH_parameters"]["config_file"])

    return config_dict


if __name__ == "__main__":

    # default parameters
    output_dir = "flow123d_sim"
    N = 2  # default number of sampling processes
    oversubscribe = False  # if there are not enough slots available
    visualize = False  # True = only visualization

    # read parameters
    len_argv = len(sys.argv)
    assert len_argv > 1, "Specify configuration yaml file!"
    if len_argv > 1:
        output_dir = sys.argv[1]
    if len_argv > 2:
        N = int(sys.argv[2])  # number of MH/DAMH chains
    if len_argv > 3:
        oversubscribe = sys.argv[3] == "oversubscribe"
        visualize = sys.argv[3] == "visualize"

    # setup paths and directories
    config_dict = setup(output_dir)
    problem_path = config_dict["bayes_config_file"]

    # run sampling
    # paths are relative to repository dir
    # paths passed to surrDAMH are absolute
    command = None
    if visualize:
        # os.error("Visualization not implemented.")
        os.chdir(config_dict["script_dir"])
        C = Configuration(N, problem_path)
        args = [str(N), problem_path, output_dir]
        if os.path.exists("surrDAMH/examples/visualization/" + C.problem_name + ".py"):
            command = "python3 surrDAMH/examples/visualization/" + C.problem_name + ".py " + " ".join(args)
        else:
            command = "python3 surrDAMH/examples/visualization/general_visualization.py " + " ".join(args)
        # local command call
        print(command)
        os.system(command)
        exit(0)
    else:
        if oversubscribe:
            opt = " --oversubscribe "
        else:
            opt = " "
        sampler = "python3 -m mpi4py surrDAMH/surrDAMH/process_SAMPLER.py " + output_dir + " "
        solver = "python3 -m mpi4py surrDAMH/surrDAMH/process_SOLVER.py " + problem_path + " " + output_dir + " "
        collector = "python3 -m mpi4py surrDAMH/surrDAMH/process_COLLECTOR.py "

        # prepare running command for local run
        # or prepare PBS script for running on Metacentrum
        if not config_dict["run_on_metacentrum"]:
            # possibly change to particular mpirun for testing
            # mpirun = "/usr/local/mpich_3.4.2/bin/mpirun"
            # mpirun = "mpirun"
            # mpirun = "mpiexec -envnone"
            mpirun = "mpiexec"
            command = mpirun + " -n " + str(N) + opt + sampler \
                      + ":" + " -n 1" + opt + solver + ":" + " -n 1" + opt + collector
        else:
            met = config_dict["metacentrum"]
            common_lines = [
                'set -x',
                '\n# absolute path to output_dir',
                'output_dir="' + output_dir + '"',
                '\n# run from the repository directory',
                'cd "' + config_dict["script_dir"] + '"',
                '\n# command for running correct docker image',
                'image=$(./endorse_fterm image)',
                'sing_command="singularity exec -B ' + met["workspace_rel"] + '/:/' + met["workspace_rel"]
                        + ' docker://$image"',
                '\n# auxiliary command for opening Python environment inside docker image',
                'bash_py="bash -c \'source ./venv/bin/activate &&"',
            ]

            # prepare PBS script
            lines = [
                '#!/bin/bash',
                '#PBS -S /bin/bash',
                '#PBS -l select=1:ncpus=' + str(N+2) + ':mem=' + met["memory"],
                '#PBS -l walltime=' + met["walltime"],
                '#PBS -q ' + met["queue"],
                '#PBS -N ' + met["name"],
                '#PBS -j oe',
                '\n',
                *common_lines,
                '\n# define surrDaMH processes',
                'sampler="$bash_py ' + sampler + '\'"',
                'solver="$bash_py ' + solver + '\'"',
                'collector="$bash_py ' + collector + '\'"',
                # '\n',
                # '# load correct MPI lib',
                # 'module load mpich-3.0.2-gcc',
                # 'which mpirun',
                # 'mpirun --version',
                '\n# get hostfile and pass into container mpiexec',
                'local_host_file="$output_dir/hostfile_$PBS_JOBID"',
                'cp $PBS_NODEFILE $local_host_file',
                '\n# finally gather the full command',
                'command="$sing_command mpiexec '
                        + '-f $local_host_file '
                        + '-launcher-exec \'$sing_command mpiexec\' '
                        + '-envnone '
                        + '-n ' + str(N) + ' $sampler : '
                        + '-n 1 $solver : '
                        + '-n 1 $collector"',
                '\n', 'echo $command', 'eval $command'
            ]
            with open("pbs_job.sh", 'w') as f:
                f.write('\n'.join(lines))

    preprocess(config_dict)

    # final command call
    if not config_dict["run_on_metacentrum"]:
        # local command call
        os.chdir(config_dict["script_dir"])
        print(command)
        os.system(command)
    else:
        # PBS script
        # os.system("qsub " + os.path.join(config_dict["work_dir"], "pbs_job.sh"))
        exit(0)
