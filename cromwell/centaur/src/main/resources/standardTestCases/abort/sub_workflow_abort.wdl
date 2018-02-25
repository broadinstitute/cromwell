import "sub_workflow_aborted_import.wdl" as sub

task let_me_run {
    command {
        echo "I'm alive !!"
    }
    runtime {
       docker: "ubuntu:latest"
    }
    output {
        Boolean done = true
    }
}

task lost_in_space {
    Boolean go
    command {
        echo "Forgotten task"
    }
    runtime {
       docker: "ubuntu:latest"
    }
}

workflow sub_workflow_abort {
    call let_me_run
    call sub.inner_abort { input: go = let_me_run.done }
    call lost_in_space { input: go = inner_abort.done }
}
