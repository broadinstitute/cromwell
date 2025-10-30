version development-1.1

task gpuTask {
  command {
    echo "No GPU :("
  }
  output {
    String msg = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
    gpu: true
  }
}

workflow wf_gpu_required_not_requested {
  call gpuTask
  output {
     String out = gpuTask.msg
  }
}
