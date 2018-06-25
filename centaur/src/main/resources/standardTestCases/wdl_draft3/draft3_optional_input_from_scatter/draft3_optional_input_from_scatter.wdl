version 1.0

workflow draft3_optional_input_from_scatter {
    scatter (x in [0]) {
        scatter (y in [0]) {
            scatter (z in [0]) {
                call input_default_not_used
            }
        }
    }

    output {
      String unpacked_out = input_default_not_used.out[0][0][0]
    }
}

task input_default_not_used {
    parameter_meta {
        greeting1: "Input with default; replaced in inputs.json"
        greeting1: "Unsupplied input, no default; set on the command line"
        greeting3: "Input with upstream default; replaced in inputs.json"
        greeting4: "Input with default; default is not replaced"
    }
    input {
        String greeting1 = "replace me"
        String? greeting2
        String greeting3 = greeting1
        String greeting4 = "hello4"
    }
    command { echo ~{greeting1} ~{default="hello2" greeting2} ~{greeting3} ~{greeting4} }
    runtime { docker: "ubuntu:latest" }
    output {
        String out = read_string(stdout())
    }
}
