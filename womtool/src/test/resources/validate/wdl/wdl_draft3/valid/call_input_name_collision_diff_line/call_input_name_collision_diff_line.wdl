version 1.0

workflow x {
    call cram
    call y as shouldntBeProblematic { input:
        cram = "asdf",
        bam = cram.scram
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
        String bam
    }
    command {
        echo "."
    }
}
