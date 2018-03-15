package wdl.draft3.transforms.wdlom2wom

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.CommandPartElement
import wdl.model.draft3.elements.CommandPartElement.StringCommandPartElement
import wom.callable.RuntimeEnvironment
import wom.expression.IoFunctionSet
import wom.{CommandPart, InstantiatedCommand}
import wom.graph.LocalName
import wom.values.WomValue

object CommandPartElementToWomCommandPart {
  def convert(commandPart: CommandPartElement): ErrorOr[CommandPart] = commandPart match {
    case s: StringCommandPartElement => WdlomWomStringCommandPart(s).valid
    case _ => s"Placeholder command part not yet implemented.".invalidNel
  }
}

case class WdlomWomStringCommandPart(stringCommandPartElement: StringCommandPartElement) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: WomValue => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] = List(InstantiatedCommand(stringCommandPartElement.value)).validNel
}
// TODO case class for Placeholder Command Part
