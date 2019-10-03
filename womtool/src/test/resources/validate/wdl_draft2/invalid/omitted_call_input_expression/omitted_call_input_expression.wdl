workflow BarcodeHg38{
    Array[File] crams

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
    File intervals
    File cram
    File ref


    command <<<
        # ... do something smart
    >>>

    runtime {

    }

    output {
        File cramFoo = "out.txt"
    }
}
