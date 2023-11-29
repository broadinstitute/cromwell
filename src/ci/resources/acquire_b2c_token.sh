#!/bin/bash

# This script is used to acquire a b2c token and save it to a temporary file, which will be
# used by the dsde-toolbox/renderCiResources to render configs.

# User must provide an output file path. 
if [ -z "$1" ]; then
    echo "Error: Must specify an output file path for the environment file."
    exit 1
fi

# Acquire a b2c token using gcloud auth
echo "Using local gcloud auth to acquire a b2c token..."
B2C_TOKEN=$(gcloud auth print-access-token)
if [ $? -eq 0 ]; then
    echo "Acquired b2c token: ${B2C_TOKEN:0:4}****"
else
    echo "Failed to acquire b2c token. Is your local shell logged into gcloud?"
fi

# Create the output directory if it doesn't exist
output_dir=$(dirname "$1")
mkdir -p "$output_dir"

# Write the token to the output file
echo "export CROMWELL_B2C_TOKEN=${B2C_TOKEN}" > "$1"