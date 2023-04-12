workflow myWorkflow {
    call myTask
}

task myTask {
    command {
        echo "hello world"
    }

    runtime {
      docker: "debian:latest"
      bootDiskSizeGb: 100
      disks: "/var/log 50 SSD"
    }
    output {
        String out = read_string(stdout())
    }
}
