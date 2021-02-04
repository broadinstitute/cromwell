task hello {
  command {
    echo "Hello world"
  }
}

workflow invalid_use_reference_disks_specification {
  call hello
}
