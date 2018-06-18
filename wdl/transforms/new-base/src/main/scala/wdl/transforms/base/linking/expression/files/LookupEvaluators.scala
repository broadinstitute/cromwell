package wdl.transforms.base.linking.expression.files

import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.FileEvaluator
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wdl.transforms.base.linking.expression.values._
import wom.expression.IoFunctionSet
import wom.types.WomType
import wom.values.{WomFile,WomValue}

object LookupEvaluators {

  implicit val identifierLookupEvaluator: FileEvaluator[IdentifierLookup] = new FileEvaluator[IdentifierLookup] {
    override def predictFilesNeededToEvaluate(a: IdentifierLookup,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val expressionMemberAccessEvaluator: FileEvaluator[ExpressionMemberAccess] = new FileEvaluator[ExpressionMemberAccess] {
    override def predictFilesNeededToEvaluate(a: ExpressionMemberAccess,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType): ErrorOr[Set[WomFile]] = {
      a.expression.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
    }
  }

  implicit val identifierMemberAccessEvaluator: FileEvaluator[IdentifierMemberAccess] = new FileEvaluator[IdentifierMemberAccess] {
    override def predictFilesNeededToEvaluate(a: IdentifierMemberAccess,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val indexAccessFileEvaluator: FileEvaluator[IndexAccess] = (a, inputs, ioFunctionSet, coerceTo) => {
    (a.expressionElement.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
    a.index.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)) mapN { _ ++ _ }
  }
}
