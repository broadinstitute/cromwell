version 1.1

# Is there a better way to test this?
import "https://raw.githubusercontent.com/broadinstitute/cromwell/fc16ebb2f3d7bc59c8ea7a49460a49cf21031600/womtool/src/test/resources/validate/biscayne/valid/relative_imports/sub_wfs/foo.wdl"

workflow http_relative_imports {
  call foo.foo_wf

  output {
    Int result = foo_wf.unpacked
  }
}
