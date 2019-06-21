version 1.0






workflow Test {

  call Echo as echo

}

task Echo {

  runtime { docker: "ubuntu:latest" }

  command {
    echo "I am prefixed with spaces"
  }

}
