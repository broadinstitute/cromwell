version 1.0

task checkMemory {
  input {
    Int mem_gb
  }
  command {
    df -h /dev/shm
  }
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
  call checkMemory as big { input: mem_gb=10 }
  output {
    String smallGb = small.out
    String bigGb = big.out
  }
}