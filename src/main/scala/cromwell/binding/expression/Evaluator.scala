package cromwell.binding.expression

import cromwell.parser.WdlParser.AstNode

import scala.language.postfixOps
import scala.util.Try

trait Evaluator {
  type T
  type LookupFunction = String => T
  type Functions = WdlFunctions[T]
  def lookup: LookupFunction
  def functions: Functions
  def evaluate(ast: AstNode): Try[T]
}

