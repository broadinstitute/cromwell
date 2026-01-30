task use_this_name {

  File f

  command {}

  meta {
    email: "skroob@spaceballs.gov"
    author: "President Skroob"
    description: "Spaceballs: The Unit Test"
  }

  output {
    File f2 = "a"
  }

  runtime {
    docker: "docker image"
  }
}
