version 1.0

import "gcpbatch_google_labels.wdl"

workflow google_labels_wrapper {
  call gcpbatch_google_labels.google_labels { input:
    expected_kvps = [ ("cromwell-sub-workflow-name", "google_labels"), ("wdl-task-name", "check_labels") ],
    check_aliases = false
  }
}
