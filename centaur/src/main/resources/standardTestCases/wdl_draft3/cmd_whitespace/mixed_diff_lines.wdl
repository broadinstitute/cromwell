version 1.0

workflow Test {

	call Echo as echo

}

task Echo {

  # Deliberately different order (tabs then spaces) than same_line
	command {
		echo "I am prefixed with tabs"
    echo "I am prefixed with spaces"
	}

}
