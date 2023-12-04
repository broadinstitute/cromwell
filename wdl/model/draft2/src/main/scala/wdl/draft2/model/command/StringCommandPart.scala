package wdl.draft2.model.command

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft2.model.Declaration
import wdl.draft2.model.expression.WdlFunctions
import wom.InstantiatedCommand
import wom.values.WomValue

case class StringCommandPart(literal: String) extends WdlCommandPart {
  override def toString: String = literal

  override def instantiate(declarations: Seq[Declaration],
                           inputsMap: Map[String, WomValue],
                           functions: WdlFunctions[WomValue],
                           valueMapper: (WomValue) => WomValue
  ): ErrorOr[List[InstantiatedCommand]] = List(InstantiatedCommand(literal)).validNel
}
