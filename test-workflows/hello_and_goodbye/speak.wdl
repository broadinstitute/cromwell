version 1.0

import "./hello.wdl" as sub

workflow speak {
  input {
    String greeting_text
  }
  call sub.hello {
    input: greeting = greeting_text
  }
}