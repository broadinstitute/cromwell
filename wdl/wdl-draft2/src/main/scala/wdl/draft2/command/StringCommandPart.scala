package wdl.draft2.command

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl._
import wdl.draft2.Declaration
import wdl.draft2.expression.WdlFunctions
import wom.InstantiatedCommand
import wom.values.WomValue

case class StringCommandPart(literal: String) extends WdlCommandPart {
  override def toString: String = literal

  override def instantiate(declarations: Seq[Declaration], inputsMap: Map[String, WomValue], functions: WdlFunctions[WomValue],
                           valueMapper: (WomValue) => WomValue): ErrorOr[List[InstantiatedCommand]] = List(InstantiatedCommand(literal)).validNel
}
