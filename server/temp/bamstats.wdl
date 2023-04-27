version 1.0

workflow BamstatsWorkflow {
    input {
        File bam_input
        Int mem_gb
    }

    call Bamstats { input: bam_input=bam_input, mem_gb=mem_gb }
}

task Bamstats {
    input {
        File bam_input
        Int mem_gb
    }

    command {
        bash /usr/local/bin/bamstats ${mem_gb} ${bam_input}
    }

    output {
        File bamstats_report = "bamstats_report.zip"
    }

    runtime {
        # Pull this docker image from the `quay.io` registry
        docker: "quay.io/collaboratory/dockstore-tool-bamstats:1.25-6_1.0"
        memory: mem_gb + "GB"
    }

    meta {
        author: "Andrew Duncan"
    }
}
