import "subworkflow.wdl" as subworkflow

workflow subworkflows_in_ifs {

  if (true) {
    Array[Int] ts = range(3)
    call subworkflow.subwf as subwfTrue { input: is = ts }
  }
  if (false) {
    Array[Int] fs = range(3)
    call subworkflow.subwf as subwfFalse { input: is = fs }
  }

  output {
    Array[Int]? tjs = subwfTrue.js
    Array[Int]? fjs = subwfFalse.js
  }
}
