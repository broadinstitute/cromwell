package wdl.draft2.command

import common.validation.ErrorOr.ErrorOr
import wdl._
import wdl.draft2.Declaration
import wdl.draft2.expression.{WdlFunctions, WdlStandardLibraryFunctions}
import wom.callable.RuntimeEnvironment
import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values.WomValue
import wom.{CommandPart, InstantiatedCommand}


trait WdlCommandPart extends CommandPart {
  def instantiate(declarations: Seq[Declaration],
                  inputsMap: Map[String, WomValue],
                  functions: WdlFunctions[WomValue],
                  valueMapper: WomValue => WomValue): ErrorOr[List[InstantiatedCommand]]

  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: WomValue => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] = {
    val wdlFunctions = WdlStandardLibraryFunctions.fromIoFunctionSet(functions)
    instantiate(Seq.empty, inputsMap.map({case (localName, value) => localName.value -> value}), wdlFunctions, valueMapper)
  }
}
