task check_if_localized_as_symlink {
  File reference_file_input
  command {
     # print true if file is a symlink, otherwise print false
     if test -h ${reference_file_input}; then echo "true"; else echo "false"; fi;
  }
  output {
    Boolean is_symlink = read_boolean(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
    backend: "Papiv2-Reference-Disk-Localization"
  }
}

workflow wf_reference_disk_test {
  call check_if_localized_as_symlink
  output {
     Boolean is_input_file_a_symlink = check_if_localized_as_symlink.is_symlink
  }
}
