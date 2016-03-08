task dosthg {
    Array[File] inputs

    command {
        echo "Nop"
    }
    runtime {
        docker: "ubuntu:latest"
    }
}
workflow lots_of_inputs {
    call dosthg
}
