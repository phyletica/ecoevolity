#! /bin/bash

set -e

echo "Creating Python 3 virtual environment for compiling documentation..."

python3 -m venv pyenv-docs
source pyenv-docs/bin/activate
pip3 install --upgrade pip
pip3 install wheel
pip3 install -r python-requirements.txt
