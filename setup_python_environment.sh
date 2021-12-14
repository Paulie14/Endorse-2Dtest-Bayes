#!/bin/bash
echo "Creating python environment.."
python3 -m venv --system-site-packages venv

source venv/bin/activate
python3 --version
which python
which pip
pip install --upgrade pip
pip -V

# already installed in image: numpy
pip install matplotlib scipy ruamel.yaml
# if we want to install our own git submodules
# pip install -e bgem

#pip freeze
deactivate

