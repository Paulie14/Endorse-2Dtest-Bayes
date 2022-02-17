#!/bin/bash
set -x

# set running on metacentrum to False
sed -i '/run_on_metacentrum:/c\run_on_metacentrum: False' config.yaml

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




#command='./endorse_fterm exec bash -c \"source ./venv/bin/activate && python3 -m mpi4py run_all.py $mcmc_config $output_dir 4\"'
#echo $command
#eval $command
#./endorse_fterm exec bash -c "\"source ./venv/bin/activate ; python3 -m mpi4py run_all.py $mcmc_config $output_dir 4\""
#./endorse_fterm exec bash -c "\"which python\""
#command="source ./venv/bin/activate && python --version"
#command="which python"

#docker run --rm -it -euid=1000 -egid=1000 -etheme=light -ewho=flow -ehome=/mnt//home/paulie -v //home/paulie:/mnt//home/paulie -w //home/paulie/Workspace/Endorse-2Dtest-Bayes -v //home/paulie/Workspace/Endorse-2Dtest-Bayes://home/paulie/Workspace/Endorse-2Dtest-Bayes -v //home/paulie/Workspace://home/paulie/Workspace flow123d/geomop:master_8c1b58980 bash -c "$command"


if [ "$run" == true ]; then
  command="source ./venv/bin/activate ; python3 -m mpi4py run_all.py $problem_path $output_dir $n_chains"
  ./endorse_fterm exec "bash -c \"$command\""
fi

# visualize
if [ "$visualize" == true ]; then
  command="source ./venv/bin/activate ; python3 run_all.py $problem_path $output_dir $n_chains visualize"
  ./endorse_fterm exec "bash -c \"$command\""
fi