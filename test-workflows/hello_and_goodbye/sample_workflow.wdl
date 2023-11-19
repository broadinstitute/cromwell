version 1.0

import "./hello.wdl" as sub

task start {
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
  command {
    echo "Starting...."
  }
  output {
    String out = read_string(stdout())
  }
}

task done {
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
  command {
    echo "Done"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow main_workflow {
  call start
  call sub.hello as hello_english {
    input: greeting = "Hello"
  }
  call sub.hello as hello_spanish {
    input: greeting = "Hola"
  }
  call sub.hello as hello_french {
    input: greeting = "Bonjour"
  }
  call done
  output {
    String speak_english = hello_english.hello_out + " " + hello_english.english_goodbye_out
  }
}