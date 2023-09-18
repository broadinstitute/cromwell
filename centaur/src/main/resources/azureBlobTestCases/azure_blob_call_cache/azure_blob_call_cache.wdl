version 1.0

workflow azure_blob_call_cache {

  input {
    File file1
    File file2
    File file3
  }

  call blob_and_acr_exercise {
    input: file1 = file1, file2 = file2, file3 = file3
  }

  output {
    String s1 = read_string(file1)
  }
}


task blob_and_acr_exercise {
  input {
    File file1
    File file2
    File file3
  }

  command {
    echo "We are not doing any real work here."
  }

  output {
    String out = read_string(stdout())
  }

  runtime {
    docker: "aednichols/public-mirror:warp-tools-wx-1263"
  }
}
