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


task readBoolean {
  command {
    dd if=/dev/zero of=file_256k bs=256k count=1
  }
  output {
     Boolean out = read_boolean("file_256k")
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task readLines {
  command {
    dd if=/dev/zero of=file_256k bs=256k count=1
  }
  output {
     Array[String] out = read_lines("file_256k")
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task readTsv {
  command {
    dd if=/dev/zero of=file_256k bs=256k count=1
  }
  output {
     Array[Array[String]] out = read_tsv("file_256k")
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task readJson {
  command {
    dd if=/dev/zero of=file_256k bs=256k count=1
  }
  output {
     Object out = read_json("file_256k")
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task readObject {
  command {
    dd if=/dev/zero of=file_256k bs=256k count=1
  }
  output {
     Object out = read_object("file_256k")
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task readObjects {
  command {
    dd if=/dev/zero of=file_256k bs=256k count=1
  }
  output {
     Array[Object] out = read_objects("file_256k")
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task readMap {
  command {
    dd if=/dev/zero of=file_256k bs=256k count=1
  }
  output {
     Map[String, String] out = read_map("file_256k")
  }
  runtime {
    docker: "ubuntu:latest"
  }
} 

task readInt {
  command {
    dd if=/dev/zero of=file_256k bs=256k count=1
  }
  output {
     Int out = read_int("file_256k")
  }
  runtime {
    docker: "ubuntu:latest"
  }
} 

task readFloat {
  command {
    dd if=/dev/zero of=file_256k bs=256k count=1
  }
  output {
     Float out = read_float("file_256k")
  }
  runtime {
    docker: "ubuntu:latest"
  }
} 

workflow read_file_limits {
  call readString
  call readBoolean
  call readLines
  call readTsv
  call readJson
  call readObject
  call readObjects
  call readMap
  call readInt
  call readFloat
}
