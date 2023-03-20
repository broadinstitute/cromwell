workflow myWorkflow {
    call myTask
    call yourTask
}

task myTask {
    command {
        echo "hello world"
    }

    runtime {
      docker: "debian:latest"
      bootDiskSizeGb: 50
    }
    output {
        String out = read_string(stdout())
    }
}
task yourTask {
    command {
        echo "goodbye cruel world"
    }

    runtime {
      docker: "debian:latest"
    }
    output {
        String out = read_string(stdout())
    }
}
