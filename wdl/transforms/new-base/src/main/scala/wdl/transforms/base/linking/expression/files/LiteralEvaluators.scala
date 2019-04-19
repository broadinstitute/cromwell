package wdl.transforms.base.linking.expression.files

import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{FileEvaluator, ValueEvaluator}
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.types.{WomCompositeType, WomSingleFileType, WomType}
import wom.values.{WomFile, WomSingleFile, WomValue}


object LiteralEvaluators {
  implicit val primitiveValueEvaluator: FileEvaluator[PrimitiveLiteralExpressionElement] = new FileEvaluator[PrimitiveLiteralExpressionElement] {
    override def predictFilesNeededToEvaluate(a: PrimitiveLiteralExpressionElement,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = Set.empty[WomFile].validNel
  }

  implicit val stringLiteralEvaluator: FileEvaluator[StringLiteral] = new FileEvaluator[StringLiteral] {
    override def predictFilesNeededToEvaluate(a: StringLiteral,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = coerceTo match {
      case WomSingleFileType => Set[WomFile](WomSingleFile(a.value)).validNel
      case _ => Set.empty[WomFile].validNel
    }
  }

  implicit val objectLiteralEvaluator: FileEvaluator[ObjectLiteral] = new FileEvaluator[ObjectLiteral] {
    override def predictFilesNeededToEvaluate(a: ObjectLiteral,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = {
      def filesInObjectField(fieldAndWomTypeTuple: (String, WomType)): ErrorOr[Set[WomFile]] = {
        val (field, womType) = fieldAndWomTypeTuple
        a.elements.get(field) match {
          case Some(fieldElement) => fieldElement.predictFilesNeededToEvaluate(inputs, ioFunctionSet, womType)
          case None => s"Invalid assignment to struct. Required field $field was not specified.".invalidNel
        }
      }

      coerceTo match {
        case WomCompositeType(mapping, _) => mapping.toList.traverse(filesInObjectField).map(_.flatten.toSet)
        case _ => a.elements.values.toList.traverse(_.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)).map(_.toSet.flatten)
      }
    }
  }

  implicit val mapLiteralEvaluator: FileEvaluator[MapLiteral] = new FileEvaluator[MapLiteral] {
    override def predictFilesNeededToEvaluate(a: MapLiteral,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = {
      a.elements.toList.flatMap { case (x,y) => List(x, y) }.traverse(_.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)).map(_.toSet.flatten)
    }
  }

  implicit val arrayLiteralEvaluator: FileEvaluator[ArrayLiteral] = new FileEvaluator[ArrayLiteral] {
    override def predictFilesNeededToEvaluate(a: ArrayLiteral,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = {
      a.elements.toList.traverse(_.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)).map(_.toSet.flatten)
    }
  }

  implicit val pairLiteralEvaluator: FileEvaluator[PairLiteral] = new FileEvaluator[PairLiteral] {
    override def predictFilesNeededToEvaluate(a: PairLiteral,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = {
      (a.left.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.right.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)) mapN { _ ++ _ }
    }
  }
}

