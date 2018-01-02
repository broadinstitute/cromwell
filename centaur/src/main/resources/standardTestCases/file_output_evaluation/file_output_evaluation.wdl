workflow file_output_evaluation_wf {

    call file_output_evaluation {}
}

task file_output_evaluation {

    command {
        echo foo > bar
        echo bar > baz
    }
    runtime {
        docker: "ubuntu:zesty"
    }
    output {
        File f = read_string("baz")
    }
}
