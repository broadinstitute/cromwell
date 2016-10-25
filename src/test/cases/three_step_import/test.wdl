import "../ps.wdl"
import "../wc.wdl"
import "../cgrep.wdl"

workflow three_step {
  call ps
  call cgrep {
    input: in_file = ps.procs
  }
  call wc {
    input: in_file = ps.procs
  }
}
