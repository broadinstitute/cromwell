task createfile {
    command {
        echo "blah" > somefile.unique.txt
        echo "blah" > someotherfile.unique.txt
    }
    runtime {
        docker: "ubuntu:latest"
    }
    output {
        File out = "somefile.unique.txt"
        File out2 = "someotherfile.unique.txt"
    }
}

task globtask {
    File inputFile
    File inputFile2
    command {
        echo "blah" > outputfile1.unique.txt
        echo "blah" > outputfile2.unique.txt
    }
    runtime {
        docker: "ubuntu:latest"
    }
    output {
        Array[File] outs = glob("*.unique.txt")
    }
}

task length {
    Array[File] array
    command {
        echo "${sep=' ' array}" | wc -w
    }
    runtime {
            docker: "ubuntu:latest"
        }
    output {
        Int size = read_int(stdout())
    }
}

workflow dontglobinputs {
    call createfile
    call globtask { input: inputFile = createfile.out, inputFile2 = createfile.out2 }
    call length{ input: array = globtask.outs }
    output {
        length.size
    }
}
