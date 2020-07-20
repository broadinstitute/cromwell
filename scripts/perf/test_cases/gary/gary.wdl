task echo_files {
  File input1
  File input2
  File input3
  File input4
  File input5
  File input6
  File input7
  File input8
  File input9
  File input10
  
  Float ref_size = size(input1, "GB") + size(input2, "GB") + size(input3, "GB") + size(input4, "GB") + size(input5, "GB") + size(input6, "GB") + size(input7, "GB") + size(input8, "GB") + size(input9, "GB") + size(input10, "GB")
  
  output {
    String out = read_string(stdout())
  }

  command {
    # sleep for 1 hour
    # sleep 3600
    
    echo "result: ${ref_size}"
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: "1"
    memory: "0.1 GB"
    preemptible: 3
    disks: "local-disk 1 HDD"
    bootDiskSizeGb: 10
  }
}

workflow echo_strings {
  scatter (i in range(30)) {
      call echo_files
  }
}
