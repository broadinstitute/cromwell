import "sub.wdl" as sub

workflow wf {
  scatter (i in range(2)) {
    call sub.wf
  }
}
