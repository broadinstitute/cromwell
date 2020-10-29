version 1.0

task t {
    input { File in_file }

    # Note: This usage of size(in_file) before passing to the command is what separates this test from the
    # attempt_to_localize_bucket_as_file test.
    Float in_file_size = size(in_file)
    command <<< echo '~{in_file}' '~{in_file_size}' >>>
    runtime { docker: "ubuntu" }
    output {
        String out = read_string(stdout())
        Float out_size = size(stdout())
    }
}

workflow attempt_to_call_size_function_on_bucket {
  call t { input: in_file = "gs://hsd_salmon_index" }
  output {
      String out = t.out
      Float out_size = t.out_size
  }
}
