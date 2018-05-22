task readString {
  command {
    dd if=/dev/zero of=file_256k bs=256k count=1
  }
  output {
     String out = read_string("file_256k")
  }
  runtime {
    docker: "ubuntu:latest"
  }
}


workflow read_file_limits {
  call readString
}
