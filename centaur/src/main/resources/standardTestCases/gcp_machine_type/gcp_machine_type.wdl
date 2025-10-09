version 1.0

workflow minimal_hello_world {
  input {
    String image = "rockylinux/rockylinux:10"
    String machine_type = "n2-standard-2"
    Int preemptible = 0
  }

  call hello_world {
    input:
      image = image,
      machine_type = machine_type,
      preemptible = preemptible
  }

  output {
    String stdout = hello_world.stdout
  }
}

task hello_world {

  input {
    String image
    String machine_type
    Int preemptible
  }

  command <<<
    cat /etc/os-release
    uname -a
    cat /proc/cpuinfo
  >>>

  runtime {
    docker: image
    gcp_machine_type: machine_type
    preemptible: preemptible
  }

  meta {
    volatile: true
  }

  output {
    String stdout = read_string(stdout())
  }
}
