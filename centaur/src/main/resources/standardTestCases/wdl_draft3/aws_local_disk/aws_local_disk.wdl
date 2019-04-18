version 1.0

workflow test {

  call local_disk {
    input:
      v = "latest"
  }
}

task local_disk {

  input {
    String v
  }

  command {
    echo "If we get here, then it worked!"
  }

  runtime {
    docker: "ubuntu:" + v
    disks: "local-disk 20 SSD"
  }
}
