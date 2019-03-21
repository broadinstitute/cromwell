version 1.0

import "d.wdl" as d

workflow counterrevolutionaries {
    call d.disproportionately

    output {
        File out = disproportionately.out
    }
}

