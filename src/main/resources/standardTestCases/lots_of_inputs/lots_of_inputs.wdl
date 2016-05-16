task dosthg {
    Array[File] inputs

    command {
        echo "Nop"
    }
    output {
       String x = read_string(stdout())
    }
    runtime {
        docker: "ubuntu:latest"
    }
}
workflow lots_of_inputs {
    call dosthg
    output{
       dosthg.x
    }
}
