version 1.0

workflow optional_declarations {
  input {
    Int? wf_int
    String? wf_string
  }

  String? coerced_wf_string = wf_int
  Int? coerced_wf_int = wf_string

  call strings { input: int = 1, string = "2" }
  call strings as strings2 # No inputs specified

  output {
    Int? s1_out_int = strings.out_int
    String? s1_out_string = strings.out_string
    Int? s1_out_coerced_int = strings.out_coerced_int
    String? s1_out_coerced_string = strings.out_coerced_string

    Int? s2_out_int = strings2.out_int
    String? s2_out_string = strings2.out_string
    Int? s2_out_coerced_int = strings2.out_coerced_int
    String? s2_out_coerced_string = strings2.out_coerced_string

    Int? out_wf_int = wf_int
    String? out_wf_string = wf_string
    String? out_coerced_wf_string = coerced_wf_string
    Int? out_coerced_wf_int = coerced_wf_int
  }
}

task strings {
  input {
    Int? int
    String? string
  }

  Int? coerced_int = string
  String? coerced_string = int

  command {
    # no-op
  }

  output {
    Int? out_int = int
    String? out_string = string
    Int? out_coerced_int = coerced_int
    String? out_coerced_string = coerced_string
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
