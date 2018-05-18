task concat_task{
    Array[String] array

    command {
        echo $${sep = ' ' array} > concat
      }
      output {
        String out = read_string("concat")
      }
    runtime {
        docker: "us.gcr.io/google-containers/ubuntu-slim:0.14"
    }
}


workflow read_file_from_bucket_with_requester_pays {
    File file
    Array[String] in_array = read_lines(file)

    call concat_task{input: array = in_array}
}
