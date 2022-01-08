version 1.0

import "b.wdl" as b

workflow anesthesiologists {
    call b.bureaucratically

    output {
        File out = bureaucratically.out
    }
}

