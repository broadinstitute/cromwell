workflow oops {
  call oopsie
}

task oopsie {
  String str
  command { echo ${str} }
  runtime { docker: docker_image }
}
