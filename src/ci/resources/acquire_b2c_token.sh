#!/bin/bash

# This script acquires a b2c token using gcloud auth and writes it to a file.
# This file is used by renderCiResources to populate the CROMWELL_B2C_TOKEN in tes_application.conf.
# It should be run as a "Before Launch" task in IntelliJ so the token is refreshed regularly.

### Setting up this script for the first time in IntelliJ? 
### 1. Open the "Run/Debug Configurations" dialog for the TES repo template.
### 2. Add a new "Before Launch" task, selecting "Run External Tool". 
###    This task should run BEFORE the renderCiResources task.
### 3. Configure the new task:
###   - Name: "Acquire b2c token"
###   - Program: /full/path/to/src/ci/resources/acquire_b2c_token.sh
###   - Arguments: ./env.temp
###   - Working directory: /full/path/to/src/ci/resources
###   - Advanced Options > Uncheck "Synchronize files after execution"
### 4. Click "OK" to save the new task. Apply and close the "Run/Debug Configurations" dialog.

# User must provide an output file path as an argument.
# This can be provided by IntelliJ as a program argument in the run configuration.
if [ -z "$1" ]; then
    echo "Error: Must specify an output file path for the environment file."
    exit 1
fi

# Acquire a b2c token using gcloud auth.
# Send stderr to /dev/null so that any gcloud auth errors/update prompts don't pollute the token output.
echo "Using local gcloud auth to acquire a b2c token..."
B2C_TOKEN=$(gcloud auth print-access-token 2>/dev/null)
if [ $? -eq 0 ]; then
    echo "Acquired b2c token: ${B2C_TOKEN:0:4}****"
else
    echo "Failed to acquire b2c token. Is your local shell logged into gcloud?"
fi

# Create the output directory if it doesn't exist
output_dir=$(dirname "$1")
mkdir -p "$output_dir"

# Write the token to the output file. dsde-toolbox will call `source` on this file that we're writing, 
# so it's expected that this file contains valid bash syntax for setting an environment variable.
echo "export CROMWELL_B2C_TOKEN=${B2C_TOKEN}" > "$1"
