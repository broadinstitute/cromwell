version 1.0

workflow wdl_v1_tests {
    scatter (x in [0]) {
        scatter (y in [0]) {
            scatter (z in [0]) {
                call input_default_not_used
            }
        }
    }
}

task input_default_not_used {
    input {
        String greeting = "hello"
        String? greeting2
        String greeting3 = greeting
    }
    command { echo ~{greeting} ~{default="unset" greeting2} ~{greeting3} }
    runtime { docker: "ubuntu:latest" }
    output {
        String out = read_string(stdout())
    }
}
