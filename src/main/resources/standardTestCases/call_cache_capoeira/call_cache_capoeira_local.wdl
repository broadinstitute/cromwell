workflow call_cache_capoeira {
    # Make some files:
    call make_files { input: ready = true }
    # Call cache the made files, creating some copies:
    call make_files as make_files_cached { input: ready = make_files.done }

    # Make sure read_files is in the call cache:
    call read_files { input: ready = true, a = make_files.bananeira, b = make_files.balanca }

    # Modify that file! Oh so crafty!
    call modify_file_sfs { input: ready = read_files.done, file_path_raw = make_files.bananeira }

    # The file inputs are different., shouldn't hit the call cache.
    call read_files as read_files_after_modify { input: ready = modify_file_sfs.done, a = make_files.bananeira, b = make_files.balanca }

    # Give the last call a chance to be writen to the call caching DB
    call sleep { input: duration = 5, ready = read_files_after_modify.done }

    # Using the copied files after modification, even if the task name is different, should still call cache to read_files:
    # ** Should call cache to read_files **
    call read_copied_files { input: ready = sleep.done, a = make_files_cached.bananeira, b = make_files_cached.balanca }

    # Task with inputs in the other order. Should not call cache:
    call read_files_swapped { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }

    # Extra whitespace and rearranged task definition.
    # Should call cache to read_files_after_modify:
    call read_files_whitespace { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }

    # Different command? Don't call cache!
    call read_files_new_command { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }

    # Different output expressions? Don't call cache!
    call read_files_new_output_expressions { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }

    # Different output names? Can't call cache (yet, at least)!
    call read_files_new_output_names { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }

    # Different non-File input? Can't call cache!
    call read_files as read_files_non_file_input_switcheroo { input: ready = !sleep.done, a = make_files.bananeira, b = make_files.balanca }

    # Different input names? Don't call cache!
    call read_files_inputs_renamed { input: ready = sleep.done, c = make_files.bananeira, d = make_files.balanca }

    # Local runtime attributes:
    call read_files_different_docker { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }
    call read_files_without_docker { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }
    call read_files_different_continueOnReturnCode { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }
    call read_files_different_continueOnReturnCode_2 { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }
    call read_files_without_continueOnReturnCode { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }
    call read_files_different_failOnStderr { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }
    call read_files_without_failOnStderr { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }
    call read_files_failOnStderr_expression { input: ready = sleep.done, a = make_files.bananeira, b = make_files.balanca }
}

task make_files {
    Boolean ready
    command {
        echo bananeira > bananeira.txt
        chmod go+w bananeira.txt
        echo balanca > balanca.txt
    }
    runtime { docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950" }
    output {
        Boolean done = true
        File bananeira = "bananeira.txt"
        File balanca = "balanca.txt"
    }
}

task modify_file_sfs {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    String file_path_raw
    # Remove if this task works without it:
    # String file_path = sub(file_path_raw, "file://", "")

    command {
        echo "ginga" > ${file_path_raw}
    }
    output {
        Boolean done = true
    }
    runtime {
        # No docker, we need this to run against the same FS as Cromwell so it can modify the right file!
        # docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        backend: "Local"
    }
}

task read_files {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File a
    File b
    command { cat ${a} ${b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

task sleep {
    Int duration
    Boolean ready
    command { sleep ${duration} }
    output { Boolean done = true }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

# Identical to read_files, apart from the name...
task read_copied_files {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File a
    File b
    command { cat ${a} ${b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

task read_files_swapped {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File a
    File b
    command { cat ${b} ${a} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

task read_files_whitespace {
    # Just like read_files, but with additional whitespace!
    Boolean   ready

    File b
    File a
    command {

                cat ${a} ${b} > out
    }
    output {
        String s = read_string("out")
        Boolean done = true
    }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

task read_files_new_command {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File a
    File b
    command {
      for i in ${a} ${b}
      do
        cat $i >> out
      done
    }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

task read_files_new_output_expressions {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File a
    File b
    command { cat ${a} ${b} > out }
    output {
        Boolean done = true
        String s = read_string("out") + "_suffix" }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

task read_files_new_output_names {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File a
    File b
    command { cat ${a} ${b} > out }
    output {
        Boolean dunn = true
        String ess = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

task read_files_inputs_renamed {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File c
    File d
    command { cat ${c} ${d} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

task read_files_different_docker {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File a
    File b
    command { cat ${a} ${b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:45b23dee08af5e43a7fea6c4cf9c25ccf269ee113168c19722f87876677c5cb2"
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

task read_files_without_docker {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File a
    File b
    command { cat ${a} ${b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

task read_files_different_continueOnReturnCode {
     # 'ready' is ignored, but useful for flow-control
     Boolean ready
     File a
     File b
     command { cat ${a} ${b} > out }
     output {
         Boolean done = true
         String s = read_string("out") }
     runtime {
         docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
         continueOnReturnCode: false
         failOnStderr: false
     }
 }

 task read_files_different_continueOnReturnCode_2 {
      # 'ready' is ignored, but useful for flow-control
      Boolean ready
      File a
      File b
      command { cat ${a} ${b} > out }
      output {
          Boolean done = true
          String s = read_string("out") }
      runtime {
          docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
          continueOnReturnCode: [ 0 ]
          failOnStderr: false
      }
  }

task read_files_without_continueOnReturnCode {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File a
    File b
    command { cat ${a} ${b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        failOnStderr: false
    }
}

task read_files_different_failOnStderr {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File a
    File b
    command { cat ${a} ${b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: true
    }
}

task read_files_without_failOnStderr {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File a
    File b
    command { cat ${a} ${b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
    }
}

task read_files_failOnStderr_expression {
    # 'ready' is ignored, but useful for flow-control
    Boolean ready
    File a
    File b
    command { cat ${a} ${b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: !ready
    }
}


