version 1.0

task dont_create_file {
    command {
      echo "No file is being generated here..."
    }
    output {
        File? nope = "nope"
    }
    runtime {
       docker: "ubuntu:latest"
    }
}

workflow wdl_optional_outputs {
    call dont_create_file
}
