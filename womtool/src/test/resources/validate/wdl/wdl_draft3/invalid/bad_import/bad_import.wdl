version 1.0

import "sub_dir/bad_import.wdl" as found_but_invalid

workflow foo {
  call found_but_invalid.call_a_bad_thing
}
