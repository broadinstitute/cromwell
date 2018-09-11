version 1.0

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
