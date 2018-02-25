
task expression {
  String version
  File memory_sizer

  command {
    echo "Hello World"
  }

  output {
    String out = read_string(stdout())
  }

  runtime {
    # Not a literal WdlString
    docker: "ubuntu:" + version

    # Uses IO functions to evaluate the 'size':
    memory: ceil(size(memory_sizer)) + "GB"
  }
}

task make_2_byte_file {
  command {
    # One byte comes from the newline:
    echo "a" > ab.txt
  }
  output {
    File f = "ab.txt"
  }
  runtime { docker: "ubuntu:latest" }
}

workflow runtime_attribute_expressions {

  call make_2_byte_file
  call expression { input: version="latest", memory_sizer = make_2_byte_file.f }

  output {
     String out = expression.out
   }
}
