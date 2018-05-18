task globbing_behavior_task{
    command {
        echo "staticFile" > staticFile.txt
        echo "staticArray" > staticArray.txt
        echo "globFile" > globFile.txt
        echo "globArray" > globArray.txt
    }

    output {
        Array[File] globArray = glob("*.txt")
        String globArrayContent = read_string(globArray[0])
    }

    runtime {
        docker: "us.gcr.io/google-containers/ubuntu-slim:0.14"
    }
}


workflow globbing_behavior_with_requester_pays{
    call globbing_behavior_task

    output {
        String globArrayContent = globbing_behavior_task.globArrayContent
    }
}
