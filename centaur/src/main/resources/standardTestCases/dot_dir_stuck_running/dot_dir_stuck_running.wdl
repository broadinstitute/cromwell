task write_path {
    command {
        mkdir out
        echo "hi" > out/hello.txt
    }

    runtime {
        docker: "ubuntu"
    }

    output {
        File dotted = "./out/hello.txt"
    }
}

task read_path {
    File dotted

    command {
        cat ${dotted}
    }

    runtime {
        docker: "ubuntu"
    }

    output {
        String out = read_string(stdout())
    }
}

# https://github.com/broadinstitute/cromwell/pull/2512
workflow dot_dir_stuck_running {
    call write_path
    call read_path as read_path_that_succeeds_but_does_not_go_to_next_call {
        input: dotted = write_path.dotted
    }

    output {
        File dotted = write_path.dotted
        String out = read_path_that_succeeds_but_does_not_go_to_next_call.out
    }
}
