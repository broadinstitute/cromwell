version 1.0

task check_if_localized_as_symlink {
  input {
    File broad_reference_file_input
    File nirvana_reference_file_input
  }
  command {
    # Print true if input is a symlink, otherwise print false.

    OUT="broad_input_symlink.txt"
    if test -h ~{broad_reference_file_input}
    then
      echo "true" > $OUT
    else
      echo "false" > $OUT
    fi

    OUT="nirvana_input_symlink.txt"
    if test -h ~{nirvana_reference_file_input}
    then
      echo "true" > $OUT
    else
      echo "false" > $OUT
    fi
  }
  output {
    Boolean is_broad_input_symlink = read_boolean("broad_input_symlink.txt")
    Boolean is_nirvana_input_symlink = read_boolean("nirvana_input_symlink.txt")
  }
  runtime {
    docker: "ubuntu:latest"
    backend: "Papiv2-Reference-Disk-Localization"
  }
}

workflow wf_reference_disk_test {
  input {
    File broad_reference_file_input
    File nirvana_reference_file_input
  }
  call check_if_localized_as_symlink {
    input:
      broad_reference_file_input = broad_reference_file_input,
      nirvana_reference_file_input = nirvana_reference_file_input
  }
  output {
    Boolean is_broad_input_file_a_symlink = check_if_localized_as_symlink.is_broad_input_symlink
    Boolean is_nirvana_input_file_a_symlink = check_if_localized_as_symlink.is_nirvana_input_symlink
  }
}
