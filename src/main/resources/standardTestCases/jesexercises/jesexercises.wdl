# Create a File result by name, and by using an engine function.
# We don't want to fail on stderr.
task url_grab {
    String URL
    command <<<
        curl -# ${URL} > news
    >>>
    output {
        File news = "news"
        File progressBar = stderr()
    }
    runtime {
        docker: "tutum/curl:latest"
        failOnStderr: false
    }
}

# Use a variety of inputs and more output variations.
task fileinfo {
    Int count
    File inFile
    command <<<
        head -n ${count} ${inFile} > ${inFile}-head
        wc -l < ${inFile}
        tail -n ${count} ${inFile} > ${inFile}-tail
    >>>
    output {
        File mirror = inFile
        File fileHead = inFile + "-head"
        Int fileLength = read_int(stdout())
        String fileTail = read_string(inFile + "-tail")
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

task cat {
    File file1
    File file2
    File file3
    String file2suffix
    String file3suffix
    command <<<
        cat ${file1} ${file2} ${file3}
        cat ${file1} > ${file1}.new
        cat ${file2} > ${file2 + file2suffix}
        cat ${file3} > ${file3}${file3suffix}
    >>>
    output {
        String catted = read_string(stdout())
        String wdlCatted = read_string(file1 + ".new") + read_string(file2 + file2suffix) + read_string(file3 + file3suffix)
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

# Redirects some stdout to stderr and some stderr to stdout.
# Also checks that JES jobs can be over 2 lines.
task outputRedirecting {
    command <<<
        logger -s "should be on stdout" 2>&1
        echo "should be on stderr" >&2
    >>>
    output {
             String stdout = read_string(stdout())
             String stderr = read_string(stderr())
         }
    runtime {
        failOnStderr: false
        docker: "ubuntu:latest"
    }
}

workflow JES_Exercises {
    File externalFileRef
    call url_grab as newsgrab {
        input:  URL = "http://www.bbc.com/news"
    }
    call fileinfo as newsinfo {
        input:  count = 20,
                inFile = newsgrab.news
    }
    call cat as Jess {
        input:  file1 = externalFileRef,
                file2 = newsgrab.progressBar,
                file2suffix="hohoho",
                file3 = newsinfo.fileHead,
                file3suffix="hahaha"
    }
    call outputRedirecting as testRedirect
    output{
       testRedirect.stdout
       testRedirect.stderr
    }
}
