task create_file {
  command {
    echo "Test file" > testfile.txt
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    File testFile = "testfile.txt"
  }
}

task read_file_with_newline_added_to_the_end_of_the_path {
  File inputFile
  File fileWithModifiedPath = inputFile + "\n"
  command {
    cat ${fileWithModifiedPath}
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow gcs_path_ending_with_newline {
  call create_file
  call read_file_with_newline_added_to_the_end_of_the_path { input: inputFile = create_file.testFile }

  output {
    String out = read_file_with_newline_added_to_the_end_of_the_path.out
  }
}
