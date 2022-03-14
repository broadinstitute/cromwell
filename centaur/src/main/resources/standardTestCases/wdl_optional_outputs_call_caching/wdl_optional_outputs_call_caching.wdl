version 1.0

task do_and_do_not_output {
    command <<<
        touch do_output.txt
    >>>
    runtime {
        docker: "ubuntu"
    }
    output {
        File? do_not_output = "do_not_output.txt"
        File? do_output = "do_output.txt"
    }
}

workflow missing_optional_output {
    call do_and_do_not_output
    output {
        File? should_be_present = do_and_do_not_output.do_output
        File? should_be_null = do_and_do_not_output.do_not_output
    }
}
