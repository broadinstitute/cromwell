version 1.0

workflow oops {
  call oopsie
}

task oopsie {
  input {
    String str
  }
  command { echo ${str} }
  runtime { docker: docker_image }
}
