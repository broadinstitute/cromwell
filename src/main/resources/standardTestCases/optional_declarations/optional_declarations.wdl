task strings {
  Int? int
  String? string

  Int? coerced_int = string
  String? coerced_string = int

  command {}
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

workflow optional_declarations {
  Int? wf_int
  String? wf_string
  String? coerced_wf_string = wf_int
  Int? coerced_wf_int = wf_string
  call strings { input: int = 1, string = "2" }
  output {
    strings.*
    Int? out_wf_int = wf_int
    String? out_wf_string = wf_string
    String? out_coerced_wf_string = coerced_wf_string
    Int? out_coerced_wf_int = coerced_wf_int
  }
}
