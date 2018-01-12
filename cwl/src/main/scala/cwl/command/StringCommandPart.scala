package cwl.command

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wom._
import wom.callable.RuntimeEnvironment
import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values.WomValue
import StringCommandPart._

case class StringCommandPart(literal: String) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] =
    // TODO CWL shellquotes by default, but this shellquotes everything. Only shellquote what should be shellquoted.
    List(InstantiatedCommand(literal.shellQuote)).validNel
}

object StringCommandPart {

  implicit class ShellQuoteHelper(val s: String) extends AnyVal {
    // Escape any double quote characters in the string and surround with double quotes. The CWL spec doesn't
    // appear to say whether quotes should be of the single or double variety, but the conformance tests clearly
    // expect double quotes to allow for interpolation of environment variables.
    def shellQuote: String = '"' + s.replaceAll("\"", "\\\"") + '"'
  }
}
