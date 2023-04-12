workflow myWorkflow {
    call myTask
}

task myTask {
    command {
        echo "hello world"
        sleep 720
    }

    runtime {
      docker: "debian:latest"
      bootDiskSizeGb: 100
      disks: "/mnt/disks/newdisk 50 SSD"
    }
    output {
        String out = read_string(stdout())
    }
}
