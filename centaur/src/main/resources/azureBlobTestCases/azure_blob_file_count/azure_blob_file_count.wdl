version 1.0

workflow azure_blob_file_count {

  input {
    File file1
  }

  output {
    String s1 = read_string(file1)
  }
}
