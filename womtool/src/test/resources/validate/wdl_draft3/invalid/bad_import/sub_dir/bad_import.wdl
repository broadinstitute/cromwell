version 1.0

import "not/a/file.wdl" as not_a_file

workflow call_a_bad_thing {
  call not_a_file.not_a_task
}
