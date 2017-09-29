package wdl.command

import wdl.expression.{WdlFunctions, WdlStandardLibraryFunctions}
import wdl.values.WdlValue
import wdl._
import wom.CommandPart
import wom.expression.IoFunctionSet

trait WdlCommandPart extends CommandPart {
  def instantiate(declarations: Seq[Declaration],
                  inputsMap: Map[String, WdlValue],
                  functions: WdlFunctions[WdlValue],
                  valueMapper: WdlValue => WdlValue): String

  override def instantiate(inputsMap: Map[String, WdlValue],
                  functions: IoFunctionSet,
                  valueMapper: WdlValue => WdlValue): String = {
    val wdlFunctions = WdlStandardLibraryFunctions.fromIoFunctionSet(functions)
    instantiate(Seq.empty, inputsMap, wdlFunctions, valueMapper)
  }
}
