version 1.0

workflow drs_basename {
    input {
        File file
    }

    String workflow_level_basename = basename(file)

    call pump_up_the_base { input:
        file = file,
        pre_basenamed = workflow_level_basename,
        basednamed_at_call_site = basename(file)
    }

    output {
        Array[String] basenames = pump_up_the_base.basenames
    }
}

task pump_up_the_base {
    input {
        File file
        String pre_basenamed
        String basednamed_at_call_site
    }

    String task_level_basename = basename(file)

    command <<<
        echo ~{pre_basenamed}
        echo ~{basednamed_at_call_site}
        echo ~{task_level_basename}
        echo ~{basename(file)}
    >>>

    output {
        Array[String] basenames = read_lines(stdout())
    }

    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa"
    }

    parameter_meta {
        # We don't actually want to use this file, just play around with its basename:
        file: { localization_optional: true }
    }
}
