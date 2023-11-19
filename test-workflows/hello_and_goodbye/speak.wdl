version 1.0

import "hello.wdl" as sub

workflow speak {
  call sub.hello as english_hello {
    input: greeting = "Hello"
  }
  call sub.hello as spanish_hello {
    input: greeting = "Hola"
  }
  call sub.hello as french_hello {
    input: greeting = "Bonjour"
  } 
}