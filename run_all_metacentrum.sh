#!/bin/bash

set -x

# set running on metacentrum to True
sed -i '/run_on_metacentrum:/c\run_on_metacentrum: True' config.yaml

# number of Markov chains
n_chains=$1
# output directory
output_dir=$2

# command
run=false
visualize=false
if [ "$3" == "visualize" ]; then
  visualize=true
elif [ "$3" == "run" ]; then
  run=true
fi

rep_dir=$(pwd)

# command for running correct docker image
image=$(./endorse_fterm image)
sing_command="singularity exec -B $rep_dir:$rep_dir docker://$image"

# auxiliary command for opening Python environment inside docker image
bash_py="bash -c 'source ./venv/bin/activate &&"


# run setup, prepare PBS script (locally, single proc)
command="$sing_command $bash_py python3 -m mpi4py run_all.py $output_dir $n_chains'"
echo $command
eval $command

qsub "$output_dir/pbs_job.sh"