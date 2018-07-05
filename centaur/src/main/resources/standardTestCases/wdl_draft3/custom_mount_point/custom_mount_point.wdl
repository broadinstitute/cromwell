version 1.0

task t {
  input {
    String v
  }


  output {
    String out1 = read_string("some_file")
    String out2 = read_string("/some/mnt/some_file")
  }
  command <<<

    echo "bazqux" > some_file

    cd /some/mnt
    echo "foobar" > some_file

  >>>
  runtime {
    docker: "ubuntu:" + v
    disks: "local-disk 20 SSD, /some/mnt 20 SSD"
  }

}
workflow custom_mount_point {

  call t {
    input:
      v = "latest"
  }

  output {
    String o2 = t.out2
    String o1 = t.out1
  }
}
