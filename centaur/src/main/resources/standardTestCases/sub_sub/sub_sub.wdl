import "sub_sub_sub.wdl" as sub_sub_sub

workflow wf {
  scatter (i in range(2)) {
    call sub_sub_sub.wf
  }
}
