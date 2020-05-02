version 1.0

task square_lots {
    input {
        Int base
        String string
    }
    command {
        python <<FIN
        for row in range(~{base}):
            print('\t'.join(['~{string}'] * ~{base}))
        FIN
    }
    runtime {
        docker: "python:latest"
    }
    output {
        Array[Array[String]] outputs = read_tsv(stdout())
    }
}

workflow megadata {
    input {
        Int base
        String string
    }
    call square_lots { input:
        base = base,
        string = string
    }
    output {
        Array[Array[String]] outputs = square_lots.outputs
    }
}


