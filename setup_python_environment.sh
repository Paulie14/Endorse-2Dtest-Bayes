#!/bin/bash
echo "[Metacentrum] Make sure you called 'source load_modules.sh' before!"
echo "Creating python environment.."
python3 -m venv venv

# source ./load_modules.sh
source venv/bin/activate
python --version
which python
which pip
pip install --upgrade pip
pip -V

# pip install wheel # error then installing bih 
pip install pyyaml attrs numpy ruamel.yaml matplotlib mpi4py
pip install -e bgem

#pip freeze
deactivate

