package wdl.command

import wdl.expression.{WdlFunctions, WdlStandardLibraryFunctions}
import wdl.values.WdlValue
import wdl._
import wom.CommandPart
import wom.graph.LocalName
import wom.expression.IoFunctionSet

trait WdlCommandPart extends CommandPart {
  def instantiate(declarations: Seq[Declaration],
                  inputsMap: Map[String, WdlValue],
                  functions: WdlFunctions[WdlValue],
                  valueMapper: WdlValue => WdlValue): String

  override def instantiate(inputsMap: Map[LocalName, WdlValue],
                  functions: IoFunctionSet,
                  valueMapper: WdlValue => WdlValue): String = {
    val wdlFunctions = WdlStandardLibraryFunctions.fromIoFunctionSet(functions)
    instantiate(Seq.empty, inputsMap.map({case (localName, value) => localName.value -> value}), wdlFunctions, valueMapper)
  }
}
