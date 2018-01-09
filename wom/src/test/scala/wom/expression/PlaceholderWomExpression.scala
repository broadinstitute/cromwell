package wom.expression

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import common.validation.ErrorOr.ErrorOr
import wom.types.WomType
import wom.values.{WomFile, WomValue}

import scala.concurrent.Future
import scala.util.Try

final case class PlaceholderWomExpression(inputs: Set[String], fixedWomType: WomType) extends WomExpression {
  override def sourceString: String = "placeholder"
  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] =
    Invalid(NonEmptyList.one(s"couldn't evaluate value from inputs $inputs\tfixedWomType\t$fixedWomType\tinputValues\t$inputValues"))
  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] =
    Valid(fixedWomType)
  override def evaluateFiles(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
    Valid(Set.empty)
}

case object PlaceholderIoFunctionSet extends IoFunctionSet {
  override def readFile(path: String): Future[String] = ???
  override def writeFile(path: String, content: String): Future[WomFile] = ???
  override def stdout(params: Seq[Try[WomValue]]) = ???
  override def stderr(params: Seq[Try[WomValue]]) = ???
  override def glob(pattern: String): Future[Seq[String]] = ???
  override def size(params: Seq[Try[WomValue]]) = ???
}
