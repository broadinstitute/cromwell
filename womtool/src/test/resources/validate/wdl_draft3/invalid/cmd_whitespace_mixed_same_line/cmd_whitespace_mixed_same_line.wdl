version 1.0

workflow Test {

	call Echo as echo

}

task Echo {

  runtime { docker: "ubuntu:latest" }

  # Deliberately different order (spaces then tabs) than diff_lines
	command {
  	echo "I am prefixed with two spaces then a tab"
	}

}
