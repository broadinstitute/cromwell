task default_runtime_attributes {
    command <<<
        echo "echo 'OH NO!'" > script.sh
        echo "exit 1" >> script.sh
        chmod +x script.sh
        ./script.sh
    >>>
    output {
        String ohno = read_string(stdout())
    }
    runtime {
        docker: "ubuntu:latest"
        # continueOnReturnCode is combined from the default_runtime_attributes in the options file.
    }
}

workflow default_runtime_attributes {
    call default_runtime_attributes
}
