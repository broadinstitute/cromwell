package cromwell.binding.values

import cromwell.binding.types.WdlArrayType

case class WdlArray(wdlType: WdlArrayType, value: Seq[WdlValue]) extends WdlValue {
  // TODO: Make sure 'value' adheres to 'wdlType'
}