
task maybe_create_file {
    Boolean do_it
    command {
        if [ "${do_it}" = true ]; then
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

task more_complicated_maybe_file_tests {
    command {
        touch yes
        echo not touching no
    }
    output {
        Array[File?] array_of_maybe_files = ["yes", "no"]
        Map[String, File?] map_with_maybe_file_values = { "yes": "yes", "no": "no" }
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


workflow wdl_optional_outputs {
    call maybe_create_file as no_create_file { input: do_it = false }
    call maybe_create_file as yes_create_file { input: do_it = true }

    call more_complicated_maybe_file_tests

}
