package cromwell.binding

import cromwell.binding.types.{WdlIntegerType, WdlStringType}

object WdlImplicits {
  implicit class WdlString(val string: String) extends AnyVal {
    def toWdlValue: WdlValue = WdlValue(string, WdlStringType)
  }

  implicit class WdlInteger(val int: Int) extends AnyVal {
    def toWdlValue: WdlValue = WdlValue(int, WdlIntegerType)
  }
}