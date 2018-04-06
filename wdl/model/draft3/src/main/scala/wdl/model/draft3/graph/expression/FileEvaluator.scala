package wdl.model.draft3.graph.expression

import cats.syntax.validated._
import cats.data.Validated.Valid
import common.validation.ErrorOr.ErrorOr
import simulacrum.typeclass
import wdl.shared.model.expression.FileEvaluatorUtil
import wom.expression.IoFunctionSet
import wom.types.WomType
import wom.values.{WomFile, WomValue}
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._

import scala.language.implicitConversions

@typeclass
trait FileEvaluator[A] {
  final def evaluateFilesNeededToEvaluate(a: A,
                                          inputs: Map[String, WomValue],
                                          ioFunctionSet: IoFunctionSet,
                                          coerceTo: WomType)
                                         (implicit valueEvaluator: ValueEvaluator[A]): ErrorOr[Set[WomFile]] = {
    a.evaluateValue(inputs, ioFunctionSet) match {
      case Valid(womValue) => FileEvaluatorUtil.findFilesToDelocalize(womValue, coerceTo).toSet.validNel
      case _ => predictFilesNeededToEvaluate(a, inputs, ioFunctionSet, coerceTo)
    }
  }

  def predictFilesNeededToEvaluate(a: A,
                                   inputs: Map[String, WomValue],
                                   ioFunctionSet: IoFunctionSet,
                                   coerceTo: WomType): ErrorOr[Set[WomFile]]
}
