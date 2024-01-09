package wdl.transforms.base.linking.expression.files

import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{FileEvaluator, ValueEvaluator}
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.types.WomType
import wom.values.{WomFile, WomValue}

object LookupEvaluators {

  implicit val identifierLookupEvaluator: FileEvaluator[IdentifierLookup] = new FileEvaluator[IdentifierLookup] {
    override def predictFilesNeededToEvaluate(a: IdentifierLookup,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType
    )(implicit
      fileEvaluator: FileEvaluator[ExpressionElement],
      valueEvaluator: ValueEvaluator[ExpressionElement]
    ): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val expressionMemberAccessEvaluator: FileEvaluator[ExpressionMemberAccess] =
    new FileEvaluator[ExpressionMemberAccess] {
      override def predictFilesNeededToEvaluate(a: ExpressionMemberAccess,
                                                inputs: Map[String, WomValue],
                                                ioFunctionSet: IoFunctionSet,
                                                coerceTo: WomType
      )(implicit
        fileEvaluator: FileEvaluator[ExpressionElement],
        valueEvaluator: ValueEvaluator[ExpressionElement]
      ): ErrorOr[Set[WomFile]] =
        a.expression.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
    }

  implicit val identifierMemberAccessEvaluator: FileEvaluator[IdentifierMemberAccess] =
    new FileEvaluator[IdentifierMemberAccess] {
      override def predictFilesNeededToEvaluate(a: IdentifierMemberAccess,
                                                inputs: Map[String, WomValue],
                                                ioFunctionSet: IoFunctionSet,
                                                coerceTo: WomType
      )(implicit
        fileEvaluator: FileEvaluator[ExpressionElement],
        valueEvaluator: ValueEvaluator[ExpressionElement]
      ): ErrorOr[Set[WomFile]] =
        Set.empty[WomFile].validNel
    }

  implicit val indexAccessFileEvaluator: FileEvaluator[IndexAccess] = new FileEvaluator[IndexAccess] {
    override def predictFilesNeededToEvaluate(a: IndexAccess,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType
    )(implicit
      fileEvaluator: FileEvaluator[ExpressionElement],
      valueEvaluator: ValueEvaluator[ExpressionElement]
    ): ErrorOr[Set[WomFile]] =
      (a.expressionElement.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
       a.index.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
      ) mapN { _ ++ _ }
  }
}
