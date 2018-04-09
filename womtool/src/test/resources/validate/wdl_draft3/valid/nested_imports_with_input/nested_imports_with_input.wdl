version draft-3

import "import_layer_1.wdl"

workflow nested_imports_with_input {
  input {
    OuterStruct os
  }

  output {
    Int unpacked = os.is.i
  }
}
