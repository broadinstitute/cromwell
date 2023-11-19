version 1.0

import "./speak.wdl" as sub

task start {
  command {
    echo "Starting...."
  }
  output {
    String out = read_string(stdout())
  }
}

task done {
  command {
    echo "Done"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow main_workflow {
  call start
  call sub.speak as hello_english {
    input: greeting_text = "Hello"
  }
  call sub.speak as hello_spanish {
    input: greeting_text = "Hola"
  }
  call sub.speak as hello_french {
    input: greeting_text = "Bonjour"
  }
  call done
}