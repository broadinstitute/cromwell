package wdl.draft3.command

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.Declaration
import wdl.draft3.expression.WdlFunctions
import wom.InstantiatedCommand
import wom.values.WomValue

case class StringCommandPart(literal: String) extends WdlCommandPart {
  override def toString: String = literal

  override def instantiate(declarations: Seq[Declaration], inputsMap: Map[String, WomValue], functions: WdlFunctions[WomValue],
                           valueMapper: (WomValue) => WomValue): ErrorOr[List[InstantiatedCommand]] = List(InstantiatedCommand(literal)).validNel
}
