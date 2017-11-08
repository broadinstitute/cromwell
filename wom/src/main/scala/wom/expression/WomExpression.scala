package wom.expression

import cats.data.Validated.Valid
import common.validation.ErrorOr.ErrorOr
import wom.types.WomType
import wom.values._

import scala.concurrent.Future
import scala.util.Try

trait WomExpression {
  def sourceString: String
  def inputs: Set[String]
  def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue]
  def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType]
  def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]]
}

/**
  * It looks and acts like an expression, but it's really just a value.
  */
final case class ValueAsAnExpression(value: WomValue) extends WomExpression {
  override def sourceString: String = value.toWomString
  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = Valid(value)
  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = Valid(value.womType)
  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] = Valid(Set.empty)
  override val inputs: Set[String] = Set.empty
}

// TODO: Flesh this out (https://github.com/broadinstitute/cromwell/issues/2521)
trait IoFunctionSet {
  def readFile(path: String): Future[String]
  def writeFile(path: String, content: String): Future[WomFile]
  def stdout(params: Seq[Try[WomValue]]): Try[WomFile]
  def stderr(params: Seq[Try[WomValue]]): Try[WomFile]
  def glob(pattern: String): Seq[String]
  def size(params: Seq[Try[WomValue]]): Try[WomFloat]
}
