task create_file_without_newline {
  command {
    echo "Test file" > file_name_part1
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    File fileWithoutNewLine = "file_name_part1"
  }
}

task create_file_with_newline {
  String inputFilePath
  String pathWithNewLine = inputFilePath + "\\\nfile_name_part2"
  command {
    gsutil cp ${inputFilePath} ${pathWithNewLine}
  }
  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk"
  }
  output {
    String outPath = pathWithNewLine
  }
}

task read_file_with_newline {
  File inputFile
  command {
    cat ${inputFile}
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow gcs_file_name_with_newline_in_the_middle {
  call create_file_without_newline
  call create_file_with_newline { input: inputFilePath = create_file_without_newline.fileWithoutNewLine }
  call read_file_with_newline { input: inputFile = create_file_with_newline.outPath }

  output {
    String out = read_file_with_newline.out
  }
}
