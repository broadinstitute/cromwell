version 1.0

workflow Test {

  call Echo as echo

}

task Echo {

  command {
echo "I am prefixed with nothing"
  }

}
