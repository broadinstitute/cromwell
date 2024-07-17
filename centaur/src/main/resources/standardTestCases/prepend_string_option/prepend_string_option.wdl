version 1.0

task prepend_string_option_task {
    input {
        String? out_prefix
        String? out_suffix
    }

    command {
        touch ${out_prefix}${'.' + out_suffix}.maf
    }

    output {
        File outFile = "${out_prefix}${'.' + out_suffix}.maf"
    }

    runtime {
        docker : "ubuntu:latest"
    }
}

workflow prepend_string_option {
    call prepend_string_option_task as both {input: out_prefix = "outPrefix", out_suffix = "outSuffix"}
    call prepend_string_option_task as first {input: out_prefix = "outPrefix"}
    call prepend_string_option_task as second {input: out_suffix = "outSuffix"}
    call prepend_string_option_task as neither
    output {
        String both_provided = basename(both.outFile)
        String first_provided = basename(first.outFile)
        String second_provided = basename(second.outFile)
        String neither_provided = basename(neither.outFile)
    }
}
