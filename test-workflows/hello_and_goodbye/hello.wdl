version 1.0

import "./goodbye.wdl" as bye

task sayHello {
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
  input {
    String say_greeting
  }
  command {
    echo "${say_greeting}"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow hello {
  input {
    String greeting
  }
  call sayHello {
    input: say_greeting = greeting
  }
  call bye.goodbye as english_goodbye {
   input: farewell = "Goodbye"
  }
  call bye.goodbye as spanish_goodbye {
    input: farewell = "Adios"
  }
  call bye.goodbye as french_goodbye {
    input: farewell = "Au revoir"
  }
  output {
    String hello_out = sayHello.out
    String english_goodbye_out = english_goodbye.out
    String spanish_goodbye_out = spanish_goodbye.out
    String french_goodbye_out = french_goodbye.out
  }
}