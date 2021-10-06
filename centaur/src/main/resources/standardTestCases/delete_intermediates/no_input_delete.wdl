version 1.0

workflow no_input_delete {
    input {
        String file_created
        File f
    }

    call simple_mirror { input: f = f }
    File f_idx = f + ".idx"

    call final_task {
        input:
            simple_mirrored_file = simple_mirror.f_out,
            index_file = f_idx
    }

    output {
        Boolean done = final_task.done
    }
}


task simple_mirror {
    meta {
        volatile: true
    }
    input {
        File f
    }
    command <<< >>>
    runtime { docker: "ubuntu:latest" }
    output {
        File f_out = f
    }
}

task final_task {
    meta {
        volatile: true
    }
    input {
        File simple_mirrored_file
        File index_file
    }
    command <<< >>>
    runtime { docker: "ubuntu:latest" }
    output {
        Boolean done = true
    }
}
