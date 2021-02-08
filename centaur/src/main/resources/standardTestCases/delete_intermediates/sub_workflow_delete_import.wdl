version 1.0

task sub_workflow_task {
    input {
    }

    command <<<
        echo "I wanna be sub deleted" > delete.txt
        echo "I think I'll go for a sub walk" > keep.txt
    >>>

    runtime {
        docker: "ubuntu"
    }

    output {
        File delete = "delete.txt"
        File keep = "keep.txt"
    }
}

workflow sub_workflow_delete_import {
    call sub_workflow_task {
        input:
    }

    output {
        File delete = sub_workflow_task.delete
        File keep = sub_workflow_task.keep
    }
}
