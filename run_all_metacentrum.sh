#!/bin/bash

set -x

# output directory for the whole simulation
output_dir="flow123d_sim"
rep_dir=$(pwd)

# command for running correct docker image
image=$(./endorse_fterm image)
sing_command="singularity exec -B $rep_dir:$rep_dir docker://$image"

# auxiliary command for opening Python environment inside docker image
bash_py="bash -c 'source ./venv/bin/activate &&"

# MCMC Bayes configuration file
mcmc_config=$(realpath config_mcmc_bayes.json)

# run setup, prepare PBS script (locally, single proc)
command="$sing_command $bash_py python3 -m mpi4py run_all.py $mcmc_config $output_dir 4'"
echo $command
eval $command

# enqueue PBS task
qsub "$output_dir/pbs_job.sh"
