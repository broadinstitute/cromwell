import "../ps.wdl"
import "../wc.wdl"
import "../cgrep.wdl"

workflow three_step {
  call ps.ps
  call cgrep.cgrep {
    input: in_file = ps.procs
  }
  call wc.wc {
    input: in_file = ps.procs
  }
}
