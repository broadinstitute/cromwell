package wdl4s.wom.expression

import wdl4s.wdl.types.WdlType
import wdl4s.wdl.values.WdlValue

import scala.concurrent.Future

trait WomExpression {
  def inputs: Set[NamedExpressionInput]
  def evaluate(variableValues: ExpressionInputs, ioFunctionSet: IoFunctionSet): WdlValue
  def womType: WdlType
}

final case class PlaceholderWomExpression(womType: WdlType) extends WomExpression {
  override def inputs: Set[NamedExpressionInput] = ???
  override def evaluate(variableValues: ExpressionInputs, ioFunctionSet: IoFunctionSet): WdlValue = ???
}

final case class NamedExpressionInput(name: String, womType: WdlType)

final case class ExpressionInputs(values: Map[NamedExpressionInput, WdlValue]) extends AnyVal

trait IoFunctionSet {
  def read_file(path: String): Future[String]
  def write_file(path: String, content: String): Future[Unit]
}

case object PlaceholderIoFunctionSet extends IoFunctionSet {
  override def read_file(path: String): Future[String] = Future.successful("35")
  override def write_file(path: String, content: String): Future[Unit] = Future.successful(())
}
