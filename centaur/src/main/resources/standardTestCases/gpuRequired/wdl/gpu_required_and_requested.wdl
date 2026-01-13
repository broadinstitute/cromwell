version 1.1

task gpuTask {
  command {
    echo "I should have a GPU :)"
  }
  output {
    String msg = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
    gpu: true
    gpuCount: 1
  }
}

workflow wf_gpu_required_and_requested {
  call gpuTask
  output {
     String out = gpuTask.msg
  }
}
