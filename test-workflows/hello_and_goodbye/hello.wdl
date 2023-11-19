version 1.0

import "./goodbye.wdl" as bye

task sayHello {
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
}