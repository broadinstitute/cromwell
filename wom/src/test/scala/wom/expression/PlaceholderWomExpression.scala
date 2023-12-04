package wom.expression

import java.util.concurrent.Executors

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import common.validation.ErrorOr.ErrorOr
import wom.types.WomType
import wom.values.WomValue

import scala.concurrent.{ExecutionContext, Future}

final case class PlaceholderWomExpression(inputs: Set[String], fixedWomType: WomType) extends WomExpression {
  override def sourceString: String = "placeholder"
  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] =
    Invalid(
      NonEmptyList.one(
        s"couldn't evaluate value from inputs $inputs\tfixedWomType\t$fixedWomType\tinputValues\t$inputValues"
      )
    )
  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] =
    Valid(fixedWomType)
  override def evaluateFiles(inputValues: Map[String, WomValue],
                             ioFunctionSet: IoFunctionSet,
                             coerceTo: WomType
  ): ErrorOr[Set[FileEvaluation]] =
    Valid(Set.empty)
}

case object DefaultSizeIoFunctionSet extends EmptyIoFunctionSet {
  val DefaultFileSize = -12345L
  implicit override def ec: ExecutionContext = ExecutionContext.fromExecutor(Executors.newFixedThreadPool(1))
  override def size(path: String): Future[Long] = Future.successful(DefaultFileSize)
}
