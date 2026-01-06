#! /bin/bash

set -e

# Make sure we leave caller back from whence they called
current_dir="$(pwd)"
clean_up() {
    cd "$current_dir"
}
trap clean_up EXIT

# Get path to directory of this script
docs_dir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

conda_env_path="${docs_dir}/conda-env.yml"
conda_env_name="eco-docs"

conda_exe="conda"

if command -v micromamba >/dev/null 2>&1
then
    conda_exe="micromamba"
elif command -v mamba >/dev/null 2>&1
then
    conda_exe="mamba"
fi

echo ""
echo "Creating docs conda environment '$conda_env_name' using command:"
echo '  '"$conda_exe" env create --name "$conda_env_name" --file "$conda_env_path"
echo ""
"$conda_exe" env create --name "$conda_env_name" --file "$conda_env_path"
