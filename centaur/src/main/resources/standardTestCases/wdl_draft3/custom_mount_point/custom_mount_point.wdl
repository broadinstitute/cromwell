version 1.0

task local_and_custom {
  input {
    String v
  }


  output {
    String out1 = read_string("some_file")
    String out2 = read_string("/some/mnt/some_file")
  }
  command <<<

    echo "local disk contents" > some_file

    cd /some/mnt
    echo "custom mount contents" > some_file

  >>>
  runtime {
    docker: "ubuntu:" + v
    disks: "local-disk 20 SSD, /some/mnt 20 SSD"
  }

}

task local_only {

  input {
    String v
  }

  command {
    echo "local disk contents" > some_file
  }

  runtime {
    docker: "ubuntu:" + v
    disks: "local-disk 20 SSD"
  }

  output {
    String out = read_string("some_file")
  }
}

task custom_only {

  input {
    String v
  }

  command {
    cd /some/mnt
    echo "custom mount contents" > some_file
  }

  runtime {
    docker: "ubuntu:" + v
    disks: "/some/mnt 20 SSD"
  }

  output {
    String out = read_string("/some/mnt/some_file")
  }
}

workflow custom_mount_point {

  call local_and_custom {
    input:
      v = "latest"
  }

  call local_only {
    input:
      v = "latest"
  }

  call custom_only {
    input:
      v = "latest"
  }

  output {
    String o1 = local_and_custom.out1
    String o2 = local_and_custom.out2
    String o3 = local_only.out
    String o4 = custom_only.out
  }
}
