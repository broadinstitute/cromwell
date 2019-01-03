version 1.0

task a {
  String rld = "rld"
  input {
    String world1 = "wo${rld}"
    String world2 = "wo~{rld}"
  }
  command {
    echo ~{"hello"} ${world1}
    echo goodbye ~{world2}
    echo ~{sep=", " true="--yes" false="--no" default="foo" world1}
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String out = read_string(stdout())
  }
}

task b {
  input {
    String world = "world"
  }
  command <<<
    echo hello ${world}
    echo goodbye ~{world}
  >>>
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String out = read_string(stdout())
  }
}
