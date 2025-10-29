version 1.0

workflow CheckAwsLabelPropagation {
  input {
    String test_message = "Testing AWS label propagation"
  }

  call CheckLabels {
    input:
      message = test_message
  }

  output {
    String result = CheckLabels.output_message
    String job_info = CheckLabels.job_info
    String job_tags = CheckLabels.job_tags
  }
}

task CheckLabels {
  input {
    String message
  }

  command <<<
    set -euo pipefail

    # Print the test message
    echo "~{message}"

    # Get the instance ID from the EC2 metadata service
    INSTANCE_ID=$(cat /var/lib/cloud/data/instance-id 2>/dev/null || echo "unknown")
    echo "Instance ID: $INSTANCE_ID"

    # Get the job ID from environment
    echo "AWS Batch Job ID: ${AWS_BATCH_JOB_ID:-unknown}"

    # Try to get tags from the job if AWS CLI is available
    if command -v aws &> /dev/null; then
      if [ -n "${AWS_BATCH_JOB_ID:-}" ]; then
        echo "Checking job tags..."
        aws batch describe-jobs --jobs "$AWS_BATCH_JOB_ID" \
          --query 'jobs[0].tags' \
          --output json > job_tags.json || echo "{}" > job_tags.json

        cat job_tags.json
      else
        echo "No AWS_BATCH_JOB_ID found"
        echo "{}" > job_tags.json
      fi
    else
      echo "AWS CLI not available"
      echo "{}" > job_tags.json
    fi

    # Output for verification
    echo "Label propagation check complete" > output.txt
  >>>

  output {
    String output_message = read_string("output.txt")
    String job_info = stdout()
    String job_tags = read_json("job_tags.json")
  }

  runtime {
    docker: "ubuntu:latest"
  }
}