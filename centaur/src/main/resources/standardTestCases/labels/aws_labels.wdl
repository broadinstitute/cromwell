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

    # Try to get tags from the instance if AWS CLI is available
    if command -v aws &> /dev/null; then
      echo "Checking instance tags..."
      aws ec2 describe-tags \
        --filters "Name=resource-id,Values=$INSTANCE_ID" \
        --query 'Tags[*].[Key,Value]' \
        --output text || echo "Could not retrieve tags"

      # Get the job ID from environment
      echo "AWS Batch Job ID: ${AWS_BATCH_JOB_ID:-unknown}"

      # Try to describe the current job
      if [ -n "${AWS_BATCH_JOB_ID:-}" ]; then
        echo "Checking job tags..."
        aws batch describe-jobs --jobs "$AWS_BATCH_JOB_ID" \
          --query 'jobs[0].tags' \
          --output json || echo "Could not retrieve job tags"
      fi
    else
      echo "AWS CLI not available"
    fi

    # Output for verification
    echo "Label propagation check complete" > output.txt
    echo "Instance: $INSTANCE_ID" >> output.txt
  >>>

  output {
    String output_message = read_string("output.txt")
    String job_info = stdout()
  }

  runtime {
    docker: "ubuntu:latest"
  }
}