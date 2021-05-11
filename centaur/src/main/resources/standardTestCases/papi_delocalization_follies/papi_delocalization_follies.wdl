version 1.0

task do_it {
    command {
        echo "foo" > file.txt
    }
    output {
        # Oops that's not the right file name. Surely Cromwell's delocalization logic will notice this and issue
        # a helpful error message...
        String value = read_string("oops.txt")
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

workflow folly {
    call do_it
}