version 1.0

import "import_me.wdl" as imported

workflow call_sub {
    input {
        String myInput
    }
    call imported.input_starved {
        input:
            feedme=myInput
    }
}