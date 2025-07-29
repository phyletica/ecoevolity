#!/usr/bin/env bash

set -e

# get location of script
this_dir="$(pwd)"

echo "Loading modules specified in 'modules-to-load.sh'..."
source "${this_dir}/modules-to-load.sh" >/dev/null 2>&1 || echo "  No modules loaded"

base_dir="$(dirname "$this_dir")"

test_env_path="${this_dir}/env-test-location.sh"
if [ -e "$test_env_path" ]
then
    echo "Getting test executable path from '$test_env_path'..."
    source "$test_env_path"
else
    ecoevolity_test_exe_path="${base_dir}/test-suite-build/test/test_ecoevolity"
fi

if [ ! -e "$ecoevolity_test_exe_path" ]
then 
    echo "ERROR: test exe '$ecoevolity_test_exe_path' does not exist."
    echo "You might need to run 'prep-tests.sh' and then re-run this script."
    exit 1
fi

test_exe_dir="$(dirname "$ecoevolity_test_exe_path")"
test_log_path="${test_exe_dir}/test.log"

if [ -e "$test_log_path" ]
then
    echo "ERROR: test log '$test_log_path' already exists."
    echo "The test executable has been run previously."
    exit 1
fi

echo "Using the following test exe path:"
echo "  $ecoevolity_test_exe_path"
echo "Directing test stderr and stdout to:"
echo "  $test_log_path"

ncl_env_path="${base_dir}/dependencies/env-ncl.sh"
if [ -e "$ncl_env_path" ]
then
    echo "Sourcing '$ncl_env_path'..."
    source "$ncl_env_path"
fi

yaml_cpp_env_path="${base_dir}/dependencies/env-yaml-cpp.sh"
# if [ -e "$yaml_cpp_env_path" ]
# then
#     echo "Sourcing '$yaml_cpp_env_path'..."
#     source "$yaml_cpp_env_path"
# fi

pll_env_path="${base_dir}/dependencies/env-pll.sh"
if [ -e "$pll_env_path" ]
then
    echo "Sourcing '$pll_env_path'..."
    source "$pll_env_path"
fi

echo ""
echo "----------------------------------------------------------------------"
echo "Running the test executable..."

# Need to change directory so the tests can find data files
cd "$test_exe_dir"
./"$(basename "$ecoevolity_test_exe_path")" > "test.log" 2>&1

echo ""
echo "Tests finished"
echo "----------------------------------------------------------------------"
