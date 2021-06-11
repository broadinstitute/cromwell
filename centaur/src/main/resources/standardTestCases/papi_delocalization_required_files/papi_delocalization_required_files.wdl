version 1.0

task do_it {
    input {
        String sync
    }
    command {
        echo "foo" > file.txt
    }
    output {
        # Intentionally not the right file name. Cromwell's file evaluation should make 'oops.txt'
        # required so this job should fail.
        String value = read_string("oops.txt")
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

task check_it {
    command <<<
      # Create the expected output.
      echo "foo" > file.txt

      # Make sure the expected output was specified as "required" in the delocalization script.
      set -euo pipefail
      grep -A 1 '"/cromwell_root/file.txt"' /cromwell_root/gcs_delocalization.sh  | tail -1 | grep '"required"'
    >>>
    output {
        String value = read_string("file.txt")
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

workflow required_files {
    call check_it
    call do_it { input: sync = check_it.value }
}
