version 1.0

import "c.wdl" as c

workflow bureaucratically {
    call c.counterrevolutionaries

    output {
        File out = counterrevolutionaries.out
    }
}

