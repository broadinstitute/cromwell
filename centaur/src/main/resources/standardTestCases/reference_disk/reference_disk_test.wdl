version 1.0

task check_if_localized_with_valid_symlink {
  input {
    File broad_reference_file_input
    File nirvana_reference_file_input
    File nirvana_reference_file_metachar_input
  }
  String broad_input_valid_symlink = "broad_input_valid_symlink.txt"
  String nirvana_input_valid_symlink = "nirvana_input_valid_symlink.txt"
  String nirvana_metachar_input_valid_symlink = "nirvana_metachar_input_valid_symlink.txt"
  command <<<
      PS4='\D{+%F %T} \w $ '
      set -o nounset -o pipefail -o xtrace

      # Echo true to stdout if the argument is a symlink pointing to an extant file, otherwise echo false.
      check_if_valid_symlink() {
          local reference_input="$1"

          if [[ -h "${reference_input}" && -f $(readlink "${reference_input}") ]]; then
              echo true
          else
              echo false
          fi
      }

      check_if_valid_symlink "~{broad_reference_file_input}" > ~{broad_input_valid_symlink}
      check_if_valid_symlink "~{nirvana_reference_file_input}" > ~{nirvana_input_valid_symlink}
      check_if_valid_symlink "~{nirvana_reference_file_metachar_input}" > ~{nirvana_metachar_input_valid_symlink}

  >>>
  output {
    Boolean is_broad_input_valid_symlink = read_boolean("~{broad_input_valid_symlink}")
    Boolean is_nirvana_input_valid_symlink = read_boolean("~{nirvana_input_valid_symlink}")
    Boolean is_nirvana_metachar_input_valid_symlink = read_boolean("~{nirvana_metachar_input_valid_symlink}")
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
  call check_if_localized_with_valid_symlink {
    input:
      broad_reference_file_input = broad_reference_file_input,
      nirvana_reference_file_input = nirvana_reference_file_input,
      nirvana_reference_file_metachar_input = nirvana_reference_file_metachar_input
  }
  output {
    Boolean is_broad_input_file_a_valid_symlink = check_if_localized_with_valid_symlink.is_broad_input_valid_symlink
    Boolean is_nirvana_input_file_a_valid_symlink = check_if_localized_with_valid_symlink.is_nirvana_input_valid_symlink
    Boolean is_nirvana_metachar_input_file_a_valid_symlink = check_if_localized_with_valid_symlink.is_nirvana_metachar_input_valid_symlink
  }
}
