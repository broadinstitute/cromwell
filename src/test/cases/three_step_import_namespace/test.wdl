import "../ps.wdl" as ns1
import "../cgrep.wdl" as ns2
import "../wc.wdl" as ns3

workflow three_step {
  call ns1.ps as a1
  call ns2.cgrep as a2 {
    input: in_file=a1.procs
  }
  call ns3.wc {
    input: in_file=a1.procs
  }
}
