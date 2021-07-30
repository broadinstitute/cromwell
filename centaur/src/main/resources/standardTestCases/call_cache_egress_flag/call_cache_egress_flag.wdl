workflow call_cache_egress_flag {
  # Make some files:
  call make_files_for_egress_flag { input: ready = true }

  # Call cache the made files
  call make_files_for_egress_flag as make_files_for_egress_flag_cached { input: ready = make_files.done }
}

task make_files_for_egress_flag {
  Boolean ready
  command {
    echo egress_flag_test > egress_flag_test.txt
    echo egress_flag_test_2 > egress_flag_test_2.txt
  }
  runtime { docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950" }
  output {
    Boolean done = true
    File egress_flag_test = "egress_flag_test.txt"
    File egress_flag_test_2 = "egress_flag_test_2.txt"
  }
}