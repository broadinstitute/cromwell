task echo_file {
    File file_not_str

    command {
        echo ${file_not_str}
    }

    runtime {
        docker: "ubuntu"
    }

    output {
        String out = read_string(stdout())
    }
}

workflow bad_file_string {
    String str = "not_a_file"

    call echo_file {
        input:
            file_not_str = str
    }
}
