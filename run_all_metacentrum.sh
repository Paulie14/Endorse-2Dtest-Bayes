#!/bin/bash

set -x

# set running on metacentrum to True
sed -i '/run_on_metacentrum:/c\run_on_metacentrum: True' config.yaml

# MCMC Bayes configuration file
problem_path=$(realpath $1)
# number of Markov chains
n_chains=$2
# output directory
output_dir=$3

# command
run=false
visualize=false
if [ "$4" == "visualize" ]; then
  visualize=true
elif [ "$4" == "run" ]; then
  run=true
fi

rep_dir=$(pwd)

# command for running correct docker image
image=$(./endorse_fterm image)
sing_command="singularity exec -B $rep_dir:$rep_dir docker://$image"

# auxiliary command for opening Python environment inside docker image
bash_py="bash -c 'source ./venv/bin/activate &&"


# run setup, prepare PBS script (locally, single proc)
command="$sing_command $bash_py python3 -m mpi4py run_all.py $problem_path $output_dir $n_chains'"
echo $command
eval $command

qsub "$output_dir/pbs_job.sh"