version development

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
