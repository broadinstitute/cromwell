version 1.0

task transubstantiation {

    runtime {
        docker: "ubuntu"
    }

    command {
        echo 'foo' > output.txt
    }

    output {
        File out = "output.txt"
    }
}

workflow disproportionately {
    call transubstantiation

    output {
        File out = transubstantiation.out
    }
}

