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
    String actual_machine_type = hello_world.actual_machine_type
    String is_preemptible = hello_world.is_preemptible
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
    curl --header "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/machine-type > actual_machine_type.txt
    curl --header "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/scheduling/preemptible > is_preemptible.txt
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
    String actual_machine_type = read_string("actual_machine_type.txt")
    String is_preemptible = read_string("is_preemptible.txt")
  }
}
