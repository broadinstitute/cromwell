package cwl.command

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wom._
import wom.callable.RuntimeEnvironment
import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values._

case class StringCommandPart(literal: String) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] =
    // TODO CWL shellquotes by default, but this shellquotes everything. Only shellquote what should be shellquoted.
    List(InstantiatedCommand(literal.shellQuote)).validNel
}
