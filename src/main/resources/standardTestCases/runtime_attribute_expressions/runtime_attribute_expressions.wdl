
task expression {
  String version

  command {
    echo "Hello world!"
    sleep 2
  }

  output {
    String out = stdout()
  }

  runtime {
    # Not a literal WdlString
    docker: "ubuntu:" + version
  }
}

workflow runtime_attribute_expressions {

  call expression { input: version="latest" }

  output {
     expression.out
   }
}
