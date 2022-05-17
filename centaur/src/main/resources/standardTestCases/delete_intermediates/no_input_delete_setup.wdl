version 1.0

import "no_input_delete.wdl" as subworkflow

workflow no_input_delete_setup {

    call makeFileAndIndex { }
    # Pass the result of making the files to the subworkflow to ensure that the makeFileAndIndex task
    # is complete before the subworkflow runs.
    call subworkflow.no_input_delete as no_input_delete { input:
        file_created = makeFileAndIndex.result,
        f = "gs://cloud-cromwell-dev-self-cleaning/cromwell_execution/no_input_delete/test_execution/no_input_delete.txt"
    }

    output {
        String done = makeFileAndIndex.result
    }
}

task makeFileAndIndex {
    meta {
        volatile: true
    }
    String outputFile = 'output.txt'
    input {
    }
    command {
        set -euo pipefail # Makes sure we fail quickly if the gsutil cp fails

        gsutil cp 'gs://cloud-cromwell-dev/cromwell_execution/no_input_delete/source_files/*.*' gs://cloud-cromwell-dev-self-cleaning/cromwell_execution/no_input_delete/test_execution/
        echo $(gsutil stat gs://cloud-cromwell-dev-self-cleaning/cromwell_execution/no_input_delete/test_execution/no_input_delete.txt) > ~{outputFile}
    }
    runtime {
      docker: "gcr.io/google.com/cloudsdktool/cloud-sdk"
    }
    output {
        String result = read_string(outputFile)
    }
}
