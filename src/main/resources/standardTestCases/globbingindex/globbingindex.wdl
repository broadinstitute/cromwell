task glob_task {
    command {
        echo "staticFile" > staticFile.txt
        echo "staticArray" > staticArray.txt
        echo "globFile" > globFile.txt
        echo "globArray" > globArray.txt
    }

    output {
        Array[File] staticFiles = ["staticFile.txt"]
        Array[File] staticRef = staticFiles
        File staticFile = staticFiles[0]
        Array[File] staticArray = ["staticArray.txt"]
        File globFile = glob("globFile.txt")[0]
        Array[File] globArray = glob("globArray.txt")

        String staticFilesContent = read_string(staticFiles[0])
        String staticRefContent = read_string(staticRef[0])
        String staticFileContent = read_string(staticFile)
        String staticArrayContent = read_string(staticArray[0])
        String globFileContent = read_string(globFile)
        String globArrayContent = read_string(globArray[0])

        Int staticFilesLength = length(staticFiles)
        Int staticRefLength = length(staticRef)
        Int staticArrayLength = length(staticArray)
        Int globArrayLength = length(globArray)
    }
}

workflow globbingindex {
    call glob_task

    output {
        String staticFilesContent = glob_task.staticFilesContent
        String staticRefContent = glob_task.staticRefContent
        String staticFileContent = glob_task.staticFileContent
        String staticArrayContent = glob_task.staticArrayContent
        String globFileContent = glob_task.globFileContent
        String globArrayContent = glob_task.globArrayContent
        Int staticFilesLength = glob_task.staticFilesLength
        Int staticRefLength = glob_task.staticRefLength
        Int staticArrayLength = glob_task.staticArrayLength
        Int globArrayLength = glob_task.globArrayLength
    }
}
