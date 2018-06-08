version 1.0

workflow wdl_v1_tests {
    scatter (x in [0]) {
        scatter (y in [0]) {
            call input_default_not_used
        }
    }
}

task input_default_not_used {
    input { String greeting = "hello" }
    command { echo ~{greeting} }
    runtime { docker: "bash" }
    output {
      String out = read_string(stdout())
    }
}
