import "https://raw.githubusercontent.com/broadinstitute/centaur/591beaf8422af7c3faf51e437a91d94d13b76eba/src/main/resources/standardTestCases/aliased_subworkflows/subworkflow.wdl" as subworkflow

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
