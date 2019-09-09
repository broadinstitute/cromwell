version 1.0

workflow dedup_localizations_papi_v2 {
  call producer
  call consumer { input: first = producer.data, second = producer.data }
  call check_log { input: out_file_path = consumer.out, log_file_name = "consumer.log" }
}

task producer {
  command {
    echo "Here is some data." > data.txt
  }

  runtime {
    docker: "ubuntu:latest"
  }

  output {
    File data = "data.txt"
  }
}

task consumer {
  input {
    File first
    File second
  }

  command {
    # noop
  }

  runtime {
    docker: "ubuntu:latest"
  }

  output {
    File out = stdout()
  }
}

task check_log {
  input {
      String out_file_path
      String log_file_name
  }
  String file_log = sub(out_file_path, "/stdout$", "/" + log_file_name)
  command {
      set -euo pipefail
      gsutil cp ~{file_log} log.txt
      set +e
      grep 'Localizing input gs://cloud-cromwell-dev-self-cleaning/cromwell_execution/travis/dedup_localizations_papi_v2/' log.txt | grep -c "data.txt"
  }
  output {
      File out = stdout()
      Int num_input_localizations = read_int(stdout())
  }
  runtime { docker: "google/cloud-sdk" }
}
