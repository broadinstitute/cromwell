workflow call_cache_egress_flag {
  
  call call_cache_egress_flag_t

  output {
    String o = call_cache_egress_flag_t.out
  }

}

task call_cache_egress_flag_t {
  String input_string

  command {
    echo input_string
  }
  runtime { docker: "ubuntu:latest" }
  output {
    String out = read_string(stdout())
  }
}
