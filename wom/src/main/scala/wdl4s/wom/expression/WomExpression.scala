package wdl4s.wom.expression

import cats.data.Validated.Valid
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wdl.types.WdlType
import wdl4s.wdl.values.WdlValue

import scala.concurrent.Future

trait WomExpression {
  def inputs: Set[String]
  def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue]
  def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType]
}

final case class PlaceholderWomExpression(inputs: Set[String], fixedWomType: WdlType) extends WomExpression {
  override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = ???
  override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] = Valid(fixedWomType)
}

// TODO: Flesh this out (https://github.com/broadinstitute/cromwell/issues/2521)
trait IoFunctionSet {
  def read_file(path: String): Future[String]
  def write_file(path: String, content: String): Future[Unit]
}

case object PlaceholderIoFunctionSet extends IoFunctionSet {
  override def read_file(path: String): Future[String] = Future.successful("35")
  override def write_file(path: String, content: String): Future[Unit] = Future.successful(())
}
