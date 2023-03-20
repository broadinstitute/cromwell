workflow myWorkflow {
    call myTask
}

task myTask {
    input {
        String name
    }
    command {
        echo "hello world"
    }

    runtime {
      docker: "ubuntu:latest"
      bootDiskSizeGb: 50
    }
    output {
        String out = read_string(stdout())
    }
}