task md5 {
    File inputFile
    command {
        echo "`date`: Running checksum on ${inputFile}..."
        md5sum ${inputFile} > md5sum.txt
        echo "`date`: Checksum is complete."
    }
    output {
        File result = "md5sum.txt"
    }
    runtime {
        docker: 'ubuntu:18.04'
        preemptible: true
    }
}

workflow fileChecksum {
    File inputFile
    call md5 { input: inputFile=inputFile}
}
