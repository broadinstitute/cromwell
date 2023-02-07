version 1.0

task check_if_localized_as_symlink {
  input {
    File broad_reference_file_input
    File nirvana_reference_file_input
    File nirvana_reference_file_metachar_input
  }
  String broad_input_symlink = "broad_input_symlink.txt"
  String nirvana_input_symlink = "nirvana_input_symlink.txt"
  String nirvana_metachar_input_symlink = "nirvana_metachar_input_symlink.txt"
  command {
    # Print true if input is a symlink, otherwise print false.
    if test -h ~{broad_reference_file_input}; then echo true; else echo false; fi > ~{broad_input_symlink}
    if test -h ~{nirvana_reference_file_input}; then echo true; else echo false; fi > ~{nirvana_input_symlink}

    # Quotes added here due to the metachar in the filename.
    if test -h "~{nirvana_reference_file_metachar_input}"; then echo true; else echo false; fi > ~{nirvana_metachar_input_symlink}
  }
  output {
    Boolean is_broad_input_symlink = read_boolean("~{broad_input_symlink}")
    Boolean is_nirvana_input_symlink = read_boolean("~{nirvana_input_symlink}")
    Boolean is_nirvana_metachar_input_symlink = read_boolean("~{nirvana_metachar_input_symlink}")
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
    File nirvana_reference_file_metachar_input
  }
  call check_if_localized_as_symlink {
    input:
      broad_reference_file_input = broad_reference_file_input,
      nirvana_reference_file_input = nirvana_reference_file_input,
      nirvana_reference_file_metachar_input = nirvana_reference_file_metachar_input
  }
  output {
    Boolean is_broad_input_file_a_symlink = check_if_localized_as_symlink.is_broad_input_symlink
    Boolean is_nirvana_input_file_a_symlink = check_if_localized_as_symlink.is_nirvana_input_symlink
    Boolean is_nirvana_metachar_input_file_a_symlink = check_if_localized_as_symlink.is_nirvana_metachar_input_symlink
  }
}
