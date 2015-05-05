package cromwell.binding

import cromwell.binding.types.WdlType

case class WdlValue(value: Any, wdlType: WdlType) {
  wdlType.checkCompatible(value)
}
