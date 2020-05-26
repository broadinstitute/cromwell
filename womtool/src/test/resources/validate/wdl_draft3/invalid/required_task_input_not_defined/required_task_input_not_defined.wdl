version 1.0

workflow misses_inputs {
    input {
        Int i
    }
    call integers {input: i=i}
}

task integers {
    input {
        Int i
        Int j
    }
    command {
        echo ~{i}, ~{j}
    }
    output {String out=stdout()}
    runtime {docker: "debian:buster-slim"}
}