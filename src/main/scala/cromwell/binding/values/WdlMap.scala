package cromwell.binding.values

import cromwell.binding.types.WdlMapType

case class WdlMap(wdlType: WdlMapType, value: Map[WdlPrimitive, WdlValue]) extends WdlValue {
  // TODO: Validate that the 'value' adheres to the 'wdlType'
}