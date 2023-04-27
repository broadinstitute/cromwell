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
      disks: "/mnt/my_mnt 50 SSD, /mnt/my_mnt2 100 HDD"
    }
    output {
        String out = read_string(stdout())
    }
}
