#!/bin/bash

set -x

# set running on metacentrum to False
sed -i '/run_on_metacentrum:/c\run_on_metacentrum: False' config.yaml

# number of chains
N=2

# output directory for the whole simulation
output_dir="flow123d_sim"
rm -rf $output_dir
rm -rf saved_samples

# MCMC Bayes configuration file
mcmc_config=$(realpath config_mcmc_bayes.json)

#command='./endorse_fterm exec bash -c \"source ./venv/bin/activate && python3 -m mpi4py run_all.py $mcmc_config $output_dir 4\"'
#echo $command
#eval $command
#./endorse_fterm exec bash -c "\"source ./venv/bin/activate ; python3 -m mpi4py run_all.py $mcmc_config $output_dir 4\""
#./endorse_fterm exec bash -c "\"which python\""
#command="source ./venv/bin/activate && python --version"
#command="which python"
command="source ./venv/bin/activate ; python3 run_all.py $mcmc_config $output_dir $N oversubscribe"
./endorse_fterm exec "bash -c \"$command\""

#docker run --rm -it -euid=1000 -egid=1000 -etheme=light -ewho=flow -ehome=/mnt//home/paulie -v //home/paulie:/mnt//home/paulie -w //home/paulie/Workspace/Endorse-2Dtest-Bayes -v //home/paulie/Workspace/Endorse-2Dtest-Bayes://home/paulie/Workspace/Endorse-2Dtest-Bayes -v //home/paulie/Workspace://home/paulie/Workspace flow123d/geomop:master_8c1b58980 bash -c "$command"

#command="source ./venv/bin/activate ; python3 -m mpi4py run_all.py $mcmc_config $output_dir $N visualize"
#./endorse_fterm exec "bash -c \"$command\""