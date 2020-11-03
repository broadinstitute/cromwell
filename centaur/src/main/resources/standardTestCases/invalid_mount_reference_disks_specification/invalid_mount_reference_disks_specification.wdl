task hello {
  command {
    echo "Hello world"
  }
}

workflow invalid_mount_reference_disks_specification {
  call hello
}
