version 1.0

import "google_labels.wdl"

workflow google_labels_wrapper {
  call google_labels.google_labels { input:
    expected_kvps = [ ("cromwell-sub-workflow-name", "google_labels"), ("wdl-task-name", "check_labels") ],
    check_aliases = false
  }
}
