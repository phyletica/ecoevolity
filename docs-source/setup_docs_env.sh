#! /bin/bash

set -e

echo "Creating Python 3 virtual environment for compiling documentation..."

python3 -m venv pyenv-docs
source pyenv-docs/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install wheel
python3 -m pip install -r python-requirements.txt
