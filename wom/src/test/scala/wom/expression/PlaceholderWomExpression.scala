package wom.expression

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import common.validation.ErrorOr.ErrorOr
import wom.types.WomType
import wom.values.{WomFile, WomSingleFile, WomValue}

import scala.concurrent.{ExecutionContext, Future}
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

trait NotImplementedIoFunctionSet extends IoFunctionSet {
  val DefaultFileSize = -12345L

  override def readFile(path: String, maxBytes: Option[Int] = None, failOnOverflow: Boolean = false): Future[String] = ???
  override def writeFile(path: String, content: String): Future[WomSingleFile] = ???
  override def copyFile(pathFrom: String, targetName: String): Future[WomSingleFile] = ???
  override def stdout(params: Seq[Try[WomValue]]) = ???
  override def stderr(params: Seq[Try[WomValue]]) = ???
  override def glob(pattern: String): Future[Seq[String]] = ???
  override def listAllFilesUnderDirectory(dirPath: String): Future[Seq[String]] = ???
  override def size(path: String): Future[Long] = Future.successful(DefaultFileSize)
  override implicit def ec: ExecutionContext = ???
}

case object PlaceholderIoFunctionSet extends NotImplementedIoFunctionSet
