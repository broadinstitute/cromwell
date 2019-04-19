version 1.0

workflow long_cmd {
    call make_string
    call echo_string {
        input:
            in_string = make_string.out_string
    }

    output {
        Int num_chars = echo_string.num_chars
    }
}

task make_string {
    # https://stackoverflow.com/questions/5349718/how-can-i-repeat-a-character-in-bash/5349842#5349842
    command <<<
        printf 'chr1 %.0s' {1..20000}
    >>>

    output {
        String out_string = read_string(stdout())
    }

    runtime {
        docker: "ubuntu"
    }
}

task echo_string {
    input {
        String in_string
    }

    command {
        set -o pipefail
        echo \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} ~{in_string} \
        | wc -m
    }

    output {
        Int num_chars = read_int(stdout())
    }

    runtime {
        docker: "ubuntu"
        memory: "4GB"
    }
}
