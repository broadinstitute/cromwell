package wdl.command

import wdl.expression.{WdlFunctions, WdlStandardLibraryFunctions}
import wdl._
import wom.CommandPart
import wom.callable.RuntimeEnvironment
import wom.graph.LocalName
import wom.expression.IoFunctionSet
import wom.values.WomValue

trait WdlCommandPart extends CommandPart {
  def instantiate(declarations: Seq[Declaration],
                  inputsMap: Map[String, WomValue],
                  functions: WdlFunctions[WomValue],
                  valueMapper: WomValue => WomValue): String

  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: WomValue => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): String = {
    val wdlFunctions = WdlStandardLibraryFunctions.fromIoFunctionSet(functions)
    instantiate(Seq.empty, inputsMap.map({case (localName, value) => localName.value -> value}), wdlFunctions, valueMapper)
  }
}
