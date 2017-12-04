import "https://raw.githubusercontent.com/broadinstitute/cromwell/f2d1c3bdd5535d9f6a997eadaf136742d86adbe5/centaur/src/main/resources/standardTestCases/aliased_subworkflows/subworkflow.wdl" as subworkflow

workflow public_http_import {

  Array[Int] ts = range(3)
  call subworkflow.subwf as subwfT { input: is = ts }

  Array[Int] fs = subwfT.js
  call subworkflow.subwf as subwfF { input: is = fs }

  output {
    Array[Int] initial = ts
    Array[Int] intermediate = subwfT.js
    Array[Int] result = subwfF.js
  }
}
