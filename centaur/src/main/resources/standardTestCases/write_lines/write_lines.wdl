# Testing write_lines() within a command block

task f2a {  # Defining a task or thing to do; Is this f2a the same as in real life (second pref.?)
  File i  # File name

  command {
    cat ${i}
  }
# Is the above a concatination?

  output {
    Array[String] out = read_lines(stdout()) # Reads each line of a file as a String and returns the lines in the file as an Array[String]
  }
  runtime {docker:"ubuntu:latest"} # Allowing use of docker image to run a task in
}

task a2f {
  Array[String] strings # Defining an array of strings named strings

  command {
    cat ${write_lines(strings)}
  }
# The above is going to concatenate lines of strings to the array strings?
  output {
    File out = stdout() # What exactly is stdout()?
    String x = read_string(out) # This string is going to be the contents of out.
  }
  runtime {docker:"ubuntu:latest"} # Please explain.
}

workflow write_lines { # The actual workflow that needs to take place
  call a2f {  # Calling the second task
    input: strings=["a","b","c","d"] # The array of strings actions are taken against.
  }

  call f2a {  # Next task called
    input: i=a2f.out # The file to be used in the f2a task is the output in a2f.
  }

  call a2f as a2f_second { # Need help explaining this. Is this like a for statement? a2f_second has all the capabilities of asf? If so, what's the point?
    input: strings=f2a.out  # Fll the array strings with the read lines of stdout().
  }
  output {
    a2f_second.x  # Read the strings from the given input. Does this get put in
  }
}