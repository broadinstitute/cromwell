package wom.expression

import cats.data.Validated.Valid
import lenthall.validation.ErrorOr.ErrorOr
import wom.types.WomType
import wom.values.{WomFile, WomFloat, WomString, WomValue}

import scala.concurrent.Future
import scala.util.Try

trait WomExpression {
  def sourceString: String
  def inputs: Set[String]
  def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue]
  def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType]
  def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]]
}

final case class PlaceholderWomExpression(inputs: Set[String], fixedWomType: WomType) extends WomExpression {
  override def sourceString: String = "placeholder"
  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = Valid(WomString("42"))
  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = Valid(fixedWomType)
  override def evaluateFiles(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] = Valid(Set.empty)
}

// TODO: Flesh this out (https://github.com/broadinstitute/cromwell/issues/2521)
trait IoFunctionSet {
  def readFile(path: String): Future[String]
  def writeFile(path: String, content: String): Future[WomFile]
  def stdout(params: Seq[Try[WomValue]]): Try[WomFile]
  def stderr(params: Seq[Try[WomValue]]): Try[WomFile]
  def glob(path: String, pattern: String): Seq[String]
  def size(params: Seq[Try[WomValue]]): Try[WomFloat]
}

case object PlaceholderIoFunctionSet extends IoFunctionSet {
  override def readFile(path: String): Future[String] = ???
  override def writeFile(path: String, content: String): Future[WomFile] = ???
  override def stdout(params: Seq[Try[WomValue]]) = ???
  override def stderr(params: Seq[Try[WomValue]]) = ???
  override def glob(path: String, pattern: String) = ???
  override def size(params: Seq[Try[WomValue]]) = ???
}
