version 1.0

import "draft3_omitting_outputs.wdl" as oops

workflow importing_omitted_outputs {
  call oops.draft3_omitting_outputs
}
