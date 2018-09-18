version 1.0

#NOTE: Please do not change the spelling for 'nb_ouptputs'. This workflow is used for Call Cache test in Perf.
#      Since there are cache entries inside test db with wrong spelling, correcting the spelling
#      will fail the test as it won't be call cached anymore!
task simpleTask {
    input {
        Int seed
        Int nb_ouptputs = 10
        Array[File] input_files
    }
    command {
        for i in `seq ~{nb_ouptputs}`; do echo $i > $i.txt; done
    }
    runtime {
        docker: "ubuntu@sha256:3f119dc0737f57f704ebecac8a6d8477b0f6ca1ca0332c7ee1395ed2c6a82be7"
    }
    output {
        Array[File] outputs = glob("*.txt")
    }
}

workflow callCachingStress {
    input {
        Int scatterSize = 1000
        Array[File] input_files
    }

    scatter(i in range(scatterSize)) {
        call simpleTask { input: input_files = input_files, seed = i }
    }
}
