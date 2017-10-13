package wom.expression

import cats.data.Validated.Valid
import lenthall.validation.ErrorOr.ErrorOr
import wom.types.WdlType
import wom.values.{WdlFile, WdlFloat, WdlString, WdlValue}

import scala.concurrent.Future
import scala.util.Try

trait WomExpression {
  def sourceString: String
  def inputs: Set[String]
  def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue]
  def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType]
  def evaluateFiles(inputTypes: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]]
}

final case class PlaceholderWomExpression(inputs: Set[String], fixedWomType: WdlType) extends WomExpression {
  override def sourceString: String = "placeholder"
  override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = Valid(WdlString("42"))
  override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] = Valid(fixedWomType)
  override def evaluateFiles(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]] = Valid(Set.empty)
}

// TODO: Flesh this out (https://github.com/broadinstitute/cromwell/issues/2521)
trait IoFunctionSet {
  def readFile(path: String): Future[String]
  def writeFile(path: String, content: String): Future[WdlFile]
  def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile]
  def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile]
  def glob(path: String, pattern: String): Seq[String]
  def size(params: Seq[Try[WdlValue]]): Try[WdlFloat]
}

case object PlaceholderIoFunctionSet extends IoFunctionSet {
  override def readFile(path: String): Future[String] = ???
  override def writeFile(path: String, content: String): Future[WdlFile] = ???
  override def stdout(params: Seq[Try[WdlValue]]) = ???
  override def stderr(params: Seq[Try[WdlValue]]) = ???
  override def glob(path: String, pattern: String) = ???
  override def size(params: Seq[Try[WdlValue]]) = ???
}
