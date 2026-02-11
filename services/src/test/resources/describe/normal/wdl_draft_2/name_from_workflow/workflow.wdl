workflow use_this_name {
}

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
