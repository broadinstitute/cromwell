package wdl4s.wdl.command

import wdl4s.wdl.expression.{WdlFunctions, WdlStandardLibraryFunctions}
import wdl4s.wdl.values.WdlValue
import wdl4s.wdl._
import wdl4s.wom.CommandPart
import wdl4s.wom.expression.IoFunctionSet

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
