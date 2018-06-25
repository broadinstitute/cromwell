#  1) Mounting a SSD to a custom location and default location with custom disk size
#  2) Write a file to that mount
#  3) Changing directory within a command shouldn't break Cromwell
#  4) Use the file that was written on the mount as an output

# ChrisTM
task t {
  String v

  command {
    echo "bazqux" > some_file

    cd /some/mnt
    echo "foobar" > some_file
  }

  output {
    String out1 = read_string("some_file")
    String out2 = read_string("/some/mnt/some_file")
  }

  runtime {
    docker: "ubuntu:" + v
    disks: "local-disk 20 SSD, /some/mnt 20 SSD"
  }
}

workflow custom_mount_point {
  call t {input: v="latest"}

  output {
     String o1 = t.out1
     String o2 = t.out2
   }
}
