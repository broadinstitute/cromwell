package wdl.draft3.transforms.linking.expression.files

import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.graph.expression.{FileEvaluator, ValueEvaluator}
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wdl.draft3.transforms.linking.expression.values._
import wom.expression.IoFunctionSet
import wom.types.WomType
import wom.values.{WomFile, WomValue}


object LiteralEvaluators {
  implicit val primitiveValueEvaluator: FileEvaluator[PrimitiveLiteralExpressionElement] = new FileEvaluator[PrimitiveLiteralExpressionElement] {
    override def predictFilesNeededToEvaluate(a: PrimitiveLiteralExpressionElement,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType): ErrorOr[Set[WomFile]] = Set.empty[WomFile].validNel
  }

  implicit val stringLiteralEvaluator: FileEvaluator[StringLiteral] = new FileEvaluator[StringLiteral] {
    override def predictFilesNeededToEvaluate(a: StringLiteral,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType): ErrorOr[Set[WomFile]] = Set.empty[WomFile].validNel
  }

  implicit val objectLiteralEvaluator: FileEvaluator[ObjectLiteral] = new FileEvaluator[ObjectLiteral] {
    override def predictFilesNeededToEvaluate(a: ObjectLiteral,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType): ErrorOr[Set[WomFile]] = {
      a.elements.values.toList.traverse(_.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)).map(_.toSet.flatten)
    }
  }

  implicit val mapLiteralEvaluator: FileEvaluator[MapLiteral] = new FileEvaluator[MapLiteral] {
    override def predictFilesNeededToEvaluate(a: MapLiteral,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType): ErrorOr[Set[WomFile]] = {
      a.elements.toList.flatMap { case (x,y) => List(x, y) }.traverse(_.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)).map(_.toSet.flatten)
    }
  }

  implicit val arrayLiteralEvaluator: FileEvaluator[ArrayLiteral] = new FileEvaluator[ArrayLiteral] {
    override def predictFilesNeededToEvaluate(a: ArrayLiteral,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType): ErrorOr[Set[WomFile]] = {
      a.elements.toList.traverse(_.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)).map(_.toSet.flatten)
    }
  }

  implicit val pairLiteralEvaluator: FileEvaluator[PairLiteral] = new FileEvaluator[PairLiteral] {
    override def predictFilesNeededToEvaluate(a: PairLiteral,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType): ErrorOr[Set[WomFile]] = {
      (a.left.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.right.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)) mapN { _ ++ _ }
    }
  }
}

