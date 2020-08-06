version 1.0

task maybe_create_file {
    input {
        Boolean do_it
    }
    command {
        if [ "~{do_it}" = true ]; then
          touch maybe
        else
          echo "No file is being generated here..."
        fi
    }
    output {
        File? maybe = "maybe"
    }
    runtime {
       docker: "ubuntu:latest"
    }
}


workflow wdl_optional_outputs {
    call maybe_create_file as no_create_file { input: do_it = false }
    call maybe_create_file as yes_create_file { input: do_it = true }
}
