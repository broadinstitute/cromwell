task empty {
  command {}

  runtime {
    docker: "ubuntu"
  }

  output {
    Int notUsed = 0
  }
}

workflow no_output_delete {
  call empty
}
