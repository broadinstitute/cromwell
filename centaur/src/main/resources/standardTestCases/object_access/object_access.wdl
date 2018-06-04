workflow object_access {
  call mk_object
  call use_object { input: obj_in = mk_object.out }
  call use_field as get_a { input: int_in = mk_object.out.a }

  Object temp = mk_object.out
  call use_field as get_b { input: int_in = temp.b }

  output {
    Array[String] lines = use_object.lines
    Int int_a = get_a.int_out
    Int int_b = get_b.int_out
  }
}

task mk_object {
  command {
    echo -e "a\tb\tc"
    echo -e "1\t2\t3"
  }
  output {
    Object out = read_object(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task use_object {
  Object obj_in
  command {
    echo ${obj_in.a}
    echo ${obj_in.b}
    echo ${obj_in.c}
  }
  output {
    Array[String] lines = read_lines(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task use_field {
  Int int_in
  command {
    echo ${int_in}
  }
  output {
    Int int_out = read_int(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
