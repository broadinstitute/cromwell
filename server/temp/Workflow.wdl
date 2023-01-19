workflow myWorkflow {
    call myTask
}

task myTask {
    command {
        echo "hello world"
    }

    runtime {
      docker: "ubuntu:latest"
    }
    output {
        String out = read_string(stdout())
    }
}