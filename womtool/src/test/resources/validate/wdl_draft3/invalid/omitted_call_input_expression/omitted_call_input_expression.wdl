version 1.0

workflow BarcodeHg38{
    input {
        Array[File] crams
    }

    scatter (cram in crams){
        call CramFoo {
            input:
                intervals = ,
                cram = cram,
                ref = ref,
        }
    }

    output {
        Array [File] cramFoos = CramFoo.cramFoo
    }
}

task CramFoo {
    input {
        File intervals
        File cram
        File ref
    }

    command <<<
        # ... do something smart
    >>>

    runtime {

    }

    output {
        File cramFoo = "out.txt"
    }
}
