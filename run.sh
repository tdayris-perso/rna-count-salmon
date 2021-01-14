#!/bin/bash
set -ei

CONDA_YAML="/mnt/beegfs/pipelines/rna-count-salmon/pipeline/rna-count-salmon/envs/workflow_flamingo.yaml"

# This function only changes echo headers
# for user's sake.
function message() {
  # Define local variables
  local status=${1}         # Either INFO, CMD, ERROR or DOC
  local message="${2:-1}"   # Your message

  # Classic switch based on status
  if [ ${status} = INFO ]; then
    echo -e "\033[1;36m@INFO:\033[0m ${message}"
  elif [ ${status} = ERROR ]; then
    echo -e "\033[41m@ERROR:\033[0m ${message}"
  elif [ ${status} = DOC ]; then
    echo -e "\033[0;33m@DOC:\033[0m ${message}"
  else
    error_handling ${LINENO} 1 "Unknown message type"
  fi
}

# This function will take error messages and exit the program
function error_handling() {
  # Geathering input parameter (message, third parameter is optionnal)
  echo -ne "\n"
  local parent_lineno="$1"
  local code="$2"
  local message="${3:-1}"

  # Checking the presence or absence of message
  if [[ -n "$message" ]] ; then
    # Case message is present
    message ERROR "Error on or near line ${parent_lineno}:\n ${message}"
    message ERROR "Exiting with status ${code}"
  else
    # Case message is not present
    message ERROR "Error on or near line ${parent_lineno}"
    message ERROR "Exiting with status ${code}"
  fi

  # Exiting with given error code
  exit "${code}"
}

function help_message() {
  message DOC "Hi, I'm functionnal only at IGR's Flamingo. Do not try to run "
  message DOC "me elsewhere. I won't work."
  echo ""
  message DOC "Thanks for using me as your script for running "
  message DOC "rna-count-salmon I'm very proud to be your script today,"
  message DOC "and I hope you'll enjoy working with me."
  echo ""
  message DOC "Every time you'll see a line starting with '@', "
  message DOC "it will be because I speak."
  message DOC "In fact, I always start my speech with :"
  message DOC "'\033[0;33m@DOC:\033[0m' when i't about my functions,"
  message DOC "'\033[1;36m@INFO:\033[0m' when it's about my intentions, "
  message DOC "'\033[41m@ERROR:\033[0m', I tell you when things go wrong."
  echo ""
  message DOC "I understand very fiew things, and here they are:"
  message DOC "-h | --help        Print this help message, then exit."
  message DOC "Otherwise, run me without any arguments and I'll do magic."
  echo ""
  message DOC "A typical command line would be:"
  message DOC "bash /path/to/run.sh"
  exit 0
}

[[ $# -gt 0 ]] && help_message

CONDA='conda'
which mamba && CONDA="mamba" || echo "Mamba not found, falling back to conda."

# Loading conda
message INFO "Sourcing conda for users who did not source it before."
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate || exit error_handling "${LINENO}" 1 "Could not source conda environment."

# Install conda environment if not installed before
message INFO "Installing environment if and only if this action is needed."
$(conda info --envs | grep "rna-count-salmon" > "/dev/null" && conda compare -n rna-count-salmon "${CONDA_YAML}") &&  message INFO "Pipeline already installed! What a chance!" || ${CONDA} env create --force -f "${CONDA_YAML}"

# Check on environment variables: if env are missing
message INFO "Loading 'rna-count-salmon' environment"
conda activate rna-count-salmon || error_handling "${LINENO}" 2 "Could not activate the environment 'rna-count-salmon'."

# then installation process did not work properly
message INFO "Running pipeline if and only if it is possible"
$(export -p | grep "RNA_COUNT_LAUNCHER" --quiet) && python3 ${RNA_COUNT_LAUNCHER} flamingo || error_handling ${LINENO} 3 "Could not locate rna-count-salmon launcher at: ${RNA_COUNT_LAUNCHER}"
