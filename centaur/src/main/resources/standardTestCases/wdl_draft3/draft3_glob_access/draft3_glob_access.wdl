version 1.0

workflow draft3_glob_access{
    call access_glob_index_task

    output {
        String globArrayContent = access_glob_index_task.globArrayContent
    }
}

task access_glob_index_task {
    command <<<
        echo "staticFile" > staticFile.txt
        echo "staticArray" > staticArray.txt
        echo "globFile" > globFile.txt
        echo "globArray" > globArray.txt
    >>>

    output {
        Array[File] globArray = glob("*.txt")
        String globArrayContent = read_string(globArray[0])
    }

    runtime {
        docker: "ubuntu:latest"
    }
}
