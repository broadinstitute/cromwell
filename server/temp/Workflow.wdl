workflow myWorkflow {
    call myTask
}

task myTask {

    command {
        echo "hello world"
    }

    runtime {
      docker: "ubuntu:latest"
      bootDiskSizeGb: 50
    }

}