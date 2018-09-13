version 1.0

workflow x {
    call cram
    call y as shouldntBeProblematic { input:
        cram = cram.scram
    }
}

task cram {
    command {
        echo "."
    }
    output {
        String scram = "."
    }
}

task y {
    input {
        String cram
    }
    command {
        echo "."
    }
}
