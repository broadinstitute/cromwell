version 1.0

task do_and_do_not_output {
    command <<<
        touch yes_output.txt
    >>>
    runtime {
        docker: "ubuntu"
    }
    output {
        File? no_output = "no_output.txt"
        File? yes_output = "yes_output.txt"
    }
}

workflow missing_optional_output {
    call do_and_do_not_output
    output {
        File? should_be_null = do_and_do_not_output.no_output
        File? should_be_present = do_and_do_not_output.yes_output
    }
}
