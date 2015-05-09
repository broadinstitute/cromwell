package cromwell.binding

import java.io.File

import cromwell.binding.types.{WdlFileType, WdlIntegerType, WdlStringType}

object WdlImplicits {
  implicit class WdlString(val string: String) extends AnyVal {
    def toWdlValue: WdlValue = WdlValue(string, WdlStringType)
  }

  implicit class WdlInteger(val int: Int) extends AnyVal {
    def toWdlValue: WdlValue = WdlValue(int, WdlIntegerType)
  }

  implicit class WdlFile(val file: File) extends AnyVal {
    def toWdlValue: WdlValue = WdlValue(file, WdlFileType)
  }
}