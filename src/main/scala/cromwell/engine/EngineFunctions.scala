package cromwell.engine

import cromwell.binding.WdlFunctions
import cromwell.binding.values.{WdlString, WdlValue}

class EngineFunctions(val ctx: TaskExecutionContext) extends WdlFunctions {
  def getFunction(name: String): WdlFunction = {

    def read_int(params: Seq[WdlValue]): WdlValue = {
      // TODO: Validate parameters

      // Just assume the first parameter is a String (example only)
      val path = params.head.asInstanceOf[WdlString].value

      // For 3-step workflow, path == "stdout"

      // Use ctx along with path to construct WdlValue
      ???
    }

    def read_string(params: Seq[WdlValue]): WdlValue = ???

    name match {
      case "read_int" => read_int
      case "read_string" => read_string
    }
  }
}
