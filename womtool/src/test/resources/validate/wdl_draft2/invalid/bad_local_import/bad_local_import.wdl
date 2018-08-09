import "not/a/workflow.wdl" as some_wdl

workflow bad_local_import {
  call some_wdl.the_task
}
