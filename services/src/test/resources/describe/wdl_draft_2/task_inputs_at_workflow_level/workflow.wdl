workflow my_workflow {

    String workflow_input

    call my_task {
        input:
            passed_from_workflow = workflow_input
    }
}

task my_task {

    String passed_from_workflow
    String task_input

    command {
        echo ${passed_from_workflow}
        echo ${task_input}
    }

}