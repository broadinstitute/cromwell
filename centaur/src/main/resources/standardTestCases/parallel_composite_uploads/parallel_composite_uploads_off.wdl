version 1.0

import "parallel_composite_uploads_lib.wdl" as lib

workflow composite_status {
  call lib.make_a_fake_bam_parallel_composite_uploads_off_in_config as make_a_fake_bam
  call lib.check_composite as check { input: bam_path = make_a_fake_bam.bam }
}
