version 1.0

task size_task {
  input {
    Float sz
   }

  command {
    echo "file has size ${sz}"
  }
  output {
    File out = stdout()
  }
  runtime {
    docker: "us.gcr.io/google-containers/ubuntu-slim:0.14"
  }
}

workflow input_from_bucket_with_requester_pays {
  input {
    File file
  }
  call size_task { input: sz = size(file) }

  output {
    String out = read_string(size_task.out)
  }
}
