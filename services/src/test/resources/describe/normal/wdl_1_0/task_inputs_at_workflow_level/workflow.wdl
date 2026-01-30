version 1.0

workflow my_workflow {

    input {
        String workflow_input
    }

    call my_task {
        input:
            passed_from_workflow = workflow_input
    }
}

task my_task {

    input {
        String passed_from_workflow
        String task_input
        String? task_input_optional
        String? task_input_optional_with_default = "a"
    }

    command {
        echo ~{passed_from_workflow}
        echo ~{task_input}
        echo ~{task_input_optional}
        echo ~{task_input_optional_with_default}
    }

}
