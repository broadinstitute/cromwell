# There is no sensible way to name this workflow, so we return the empty string
task not_this_name {
  command {}

  runtime {
    docker: "docker image"
  }
}

task nor_this_one {
  command {}

  runtime {
    docker: "expression" + "based" + "image"
  }
}
