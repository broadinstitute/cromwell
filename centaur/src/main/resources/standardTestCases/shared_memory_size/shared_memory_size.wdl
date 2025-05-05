version 1.0

task checkMemory {
  input {
    Int mem_gb
  }
  command <<<
    df -h /dev/shm | tail -1 | awk '{print $2}'
  >>>
  output {
    String out = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
    memory: mem_gb + "GB"
  }
}

workflow shared_memory_size {
  call checkMemory as small { input: mem_gb=5 }
  call checkMemory as medium { input: mem_gb=10 }
  call checkMemory as big { input: mem_gb=20 }
  output {
    String smallGb = small.out
    String mediumGb = medium.out
    String bigGb = big.out
  }
}