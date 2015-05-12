package cromwell.binding

import java.io.File

import cromwell.binding.values._
import cromwell.binding.values.WdlValue

object WdlImplicits {
  implicit class WdlStringImplicit(val string: String) extends AnyVal {
    def toWdlValue: WdlValue = WdlString(string)
  }

  implicit class WdlIntegerImplicit(val int: Int) extends AnyVal {
    def toWdlValue: WdlValue = WdlInteger(int)
  }

  implicit class WdlFileImplicit(val file: File) extends AnyVal {
    def toWdlValue: WdlValue = WdlFile(file)
  }
}