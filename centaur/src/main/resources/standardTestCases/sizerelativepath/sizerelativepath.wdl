version 1.0

task print_size {
    input {
        File file
    }
    Int bytes = size(file)

    command {
        echo ~{bytes}
    }

    output {
        String out = stdout()
    }

    runtime {docker: "ubuntu:latest"}
}

workflow sizerelativepath {
    input {
        File file = "LICENSE.txt"
    }

    call print_size {
        input:
            file = file
    }
    output {
        String size_string = print_size.out
    }
}