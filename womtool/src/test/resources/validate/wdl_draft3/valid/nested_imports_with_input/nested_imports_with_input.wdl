version 1.0

import "import_layer_1.wdl" alias OuterStruct as OoterStraact

workflow nested_imports_with_input {
  input {
    OoterStraact os
  }

  output {
    Int unpacked = os.is.i
  }
}
