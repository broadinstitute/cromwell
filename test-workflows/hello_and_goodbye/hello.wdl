version 1.0

import "goodbye.wdl" as sub

task sayHello {
  input {
    String greeting
  }
  command {
    echo "${greeting}"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow hello{
  input {
    String greeting
  }
  call sayHello {
    input: greeting = greeting
  }
  call sub.goodbye as english_goodbye {
   input: farewell = "Goodbye"
  }
  call sub.goodbye as spanish_goodbye {
    input: farewell = "Adios"
  }
  call sub.goodbye as french_goodbye {
    input: farewell = "Au revoir"
  }
}