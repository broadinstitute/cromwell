version 1.0

workflow call_cache_cha_cha {
    input {
        Boolean modify_file_using_cloud
        String modify_file_command_prefix
        String modify_file_docker
    }

    # Make some files:
    call make_files { input: ready = true }
    # Call cache the made files, creating some copies:
    call make_files as make_files_cached { input: ready = make_files.done }

    # Make sure read_files is in the call cache:
    call read_files { input: ready = true, a = make_files.enganadora, b = make_files.alardoso }

    # Modify that file! Oh so crafty!
    if (modify_file_using_cloud) {
        call modify_file_cloud {
            input: ready = read_files.done,
                file_path_raw = make_files.enganadora,
                modify_file_command_prefix = modify_file_command_prefix,
                modify_file_docker = modify_file_docker
        }
    }
    if (!modify_file_using_cloud) {
        call modify_file_local {
            input: ready = read_files.done,
                file_path_raw = make_files.enganadora,
                modify_file_command_prefix = modify_file_command_prefix
        }
    }

    Boolean modify_file_done = select_first([modify_file_cloud.done, modify_file_local.done])

    # The file inputs are different., shouldn't hit the call cache.
    call read_files as read_files_after_modify {
        input: ready = modify_file_done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }

    # Give the last call a chance to be writen to the call caching DB
    call sleep as sleep_before_restart { input: duration = 5, ready = read_files_after_modify.done }

    call sleep as sleep_during_restart { input: duration = 60, ready = sleep_before_restart.done }

    call sleep as sleep_after_restart { input: duration = 5, ready = sleep_during_restart.done }

    # Using the copied files after modification, even if the task name is different, should still call cache to read_files:
    # ** Should call cache to read_files **
    call read_copied_files {
        input: ready = sleep_after_restart.done,
            a = make_files_cached.enganadora,
            b = make_files_cached.alardoso
    }

    # Task with inputs in the other order. Should not call cache:
    call read_files_swapped {
        input: ready = sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }

    # Extra whitespace and rearranged task definition.
    # Should call cache to read_files_after_modify:
    call read_files_whitespace {
        input: ready = sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }

    # Different command? Don't call cache!
    call read_files_new_command {
        input: ready = sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }

    # Different output expressions? Don't call cache!
    call read_files_new_output_expressions {
        input: ready = sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }

    # Different output names? Can't call cache (yet, at least)!
    call read_files_new_output_names {
        input: ready = sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }

    # Different non-File input? Can't call cache!
    call read_files as read_files_non_file_input_switcheroo {
        input: ready = !sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }

    # Different input names? Don't call cache!
    call read_files_inputs_renamed {
        input: ready = sleep_after_restart.done,
            c = make_files.enganadora,
            d = make_files.alardoso
    }

    # Like read_files but with 2 Array[File] instead of File
    call read_array_files {
        input: ready = sleep_after_restart.done,
            a = [make_files.enganadora],
            b = [make_files.alardoso]
    }

    # Give the last call a chance to be writen to the call caching DB
    call sleep as sleep_some_more { input: duration = 5, ready = read_array_files.done }

    # Same as read_array_files except a is empty and b contains both, in the same order. Don't call cache!
    call read_array_files as read_array_files_rearranged {
        input: ready = sleep_some_more.done,
            a = [],
            b = [make_files.enganadora,
            make_files.alardoso]
    }

    # Local runtime attributes:
    call read_files_different_docker {
        input: ready = sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }
    call read_files_different_continueOnReturnCode {
        input: ready = sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }
    call read_files_different_continueOnReturnCode_2 {
        input: ready = sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }
    call read_files_without_continueOnReturnCode {
        input: ready = sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }
    call read_files_different_failOnStderr {
        input: ready = sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }
    call read_files_without_failOnStderr {
        input: ready = sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }
    call read_files_failOnStderr_expression {
        input: ready = sleep_after_restart.done,
            a = make_files.enganadora,
            b = make_files.alardoso
    }
}

task make_files {
    input {
        Boolean ready
    }
    command {
        echo enganadora > enganadora.txt
        chmod go+w enganadora.txt
        echo alardoso > alardoso.txt
    }
    output {
        Boolean done = true
        File enganadora = "enganadora.txt"
        File alardoso = "alardoso.txt" }
    runtime { docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950" }
}

task modify_file_local {
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        String file_path_raw
        String modify_file_command_prefix
        # Remove if this task works without it:
        # String file_path = sub(file_path_raw, "file://", "")
    }
    command {
        ~{modify_file_command_prefix} ~{file_path_raw}
    }
    output { Boolean done = true }
    runtime {
        # No docker, we need this to run against the same FS as Cromwell so it can modify the right file!
        # docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        backend: "Local"
    }
}

task modify_file_cloud {
    input {
        Boolean ready
        String file_path_raw
        String modify_file_command_prefix
        String modify_file_docker
    }
    command {
        ~{modify_file_command_prefix} ~{file_path_raw}
    }
    output {
        Boolean done = true }
    runtime {
        docker: modify_file_docker
    }
}

task read_files {
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        File a
        File b
    }
    command { cat ~{a} ~{b} > out }
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
    input {
        Int duration
        Boolean ready
    }
    command { sleep ~{duration} }
    output { Boolean done = true }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

# Identical to read_files, apart from the name...
task read_copied_files {
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        File a
        File b
    }
    command { cat ~{a} ~{b} > out }
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
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        File a
        File b
    }
    command { cat ~{b} ~{a} > out }
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
    input {
        # Just like read_files, but with additional whitespace!
        Boolean   ready

        File b
        File a
    }
    command {

                cat ~{a} ~{b} > out
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

task read_files_new_command {
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        File a
        File b
    }
    command {
      for i in ~{a} ~{b}
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
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        File a
        File b
    }
    command { cat ~{a} ~{b} > out }
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
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        File a
        File b
    }
    command { cat ~{a} ~{b} > out }
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
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        File c
        File d
    }
    command { cat ~{c} ~{d} > out }
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
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        File a
        File b
    }
    command { cat ~{a} ~{b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:45b23dee08af5e43a7fea6c4cf9c25ccf269ee113168c19722f87876677c5cb2"
        continueOnReturnCode: 0
        failOnStderr: false
    }
}

task read_files_different_continueOnReturnCode {
    input {
         # 'ready' is ignored, but useful for flow-control
         Boolean ready
         File a
         File b
    }
    command { cat ~{a} ~{b} > out }
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
    input {
         # 'ready' is ignored, but useful for flow-control
         Boolean ready
         File a
         File b
    }
    command { cat ~{a} ~{b} > out }
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
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        File a
        File b
    }
    command { cat ~{a} ~{b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        failOnStderr: false
    }
}

task read_files_different_failOnStderr {
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        File a
        File b
    }
    command { cat ~{a} ~{b} > out }
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
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        File a
        File b
    }
    command { cat ~{a} ~{b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
    }
}

task read_files_failOnStderr_expression {
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        File a
        File b
    }
    command { cat ~{a} ~{b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: !ready
    }
}

task read_array_files {
    input {
        # 'ready' is ignored, but useful for flow-control
        Boolean ready
        Array[File] a
        Array[File] b
    }
    command { cat ~{sep = " " a} ~{sep = " " b} > out }
    output {
        Boolean done = true
        String s = read_string("out") }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        continueOnReturnCode: 0
        failOnStderr: !ready
    }
}
