version 1.0
task echo {
    input {
        String s
    }

    command {
        echo "~{s}"
    }

    output {
        String who = read_string(stdout())
    }
}
