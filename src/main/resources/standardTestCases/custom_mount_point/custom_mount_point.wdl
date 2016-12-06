#  1) Mounting a SSD to a custom location
#  2) Write a file to that mount
#  3) Changing directory within a command shouldn't break Cromwell
#  4) Use the file that was written on the mount as an output

task t {
  String version

  command {
    cd /some/mnt
    echo "foobar" > some_file
    sleep 2
  }

  output {
    String out2 = read_string("/some/mnt/some_file")
  }

  runtime {
    docker: "ubuntu:" + version
    disks: "local-disk 20 SSD, /some/mnt 20 SSD"
  }
}

workflow custom_mount_point {
  call t {input: version="latest"}

  output {
     t.out2
   }
}
