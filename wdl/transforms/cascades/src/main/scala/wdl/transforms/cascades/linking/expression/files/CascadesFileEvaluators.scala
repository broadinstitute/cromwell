package wdl.transforms.cascades.linking.expression.files

import cats.implicits.{catsSyntaxValidatedId, toTraverseOps}
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{
  AsMap,
  AsPairs,
  CollectByKey,
  Keys,
  Max,
  Min,
  Quote,
  Sep,
  SQuote,
  StructLiteral,
  SubPosix,
  Suffix,
  Unzip
}
import wdl.model.draft3.graph.expression.{FileEvaluator, ValueEvaluator}
import wdl.transforms.base.linking.expression.files.EngineFunctionEvaluators
import wdl.transforms.base.linking.expression.files.EngineFunctionEvaluators.{
  threeParameterFunctionPassthroughFileEvaluator,
  twoParameterFunctionPassthroughFileEvaluator
}
import wom.expression.IoFunctionSet
import wom.types.{WomCompositeType, WomType}
import wom.values.{WomFile, WomValue}

object cascadesFileEvaluators {

  implicit val keysFileEvaluator: FileEvaluator[Keys] = EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val asMapFileEvaluator: FileEvaluator[AsMap] =
    EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val asPairsFileEvaluator: FileEvaluator[AsPairs] =
    EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val collectByKeyFileEvaluator: FileEvaluator[CollectByKey] =
    EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val unzipFunctionEvaluator: FileEvaluator[Unzip] =
    EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator

  implicit val sepFunctionEvaluator: FileEvaluator[Sep] = twoParameterFunctionPassthroughFileEvaluator[Sep]
  implicit val subPosixFunctionEvaluator: FileEvaluator[SubPosix] =
    threeParameterFunctionPassthroughFileEvaluator[SubPosix]
  implicit val suffixFunctionEvaluator: FileEvaluator[Suffix] = twoParameterFunctionPassthroughFileEvaluator[Suffix]
  implicit val quoteFunctionEvaluator: FileEvaluator[Quote] =
    EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val sQuoteFunctionEvaluator: FileEvaluator[SQuote] =
    EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator

  implicit val minFunctionEvaluator: FileEvaluator[Min] = twoParameterFunctionPassthroughFileEvaluator[Min]
  implicit val maxFunctionEvaluator: FileEvaluator[Max] = twoParameterFunctionPassthroughFileEvaluator[Max]

  implicit val structLiteralEvaluator: FileEvaluator[StructLiteral] = new FileEvaluator[StructLiteral] {
    override def predictFilesNeededToEvaluate(a: StructLiteral,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType
    )(implicit
      fileEvaluator: FileEvaluator[ExpressionElement],
      valueEvaluator: ValueEvaluator[ExpressionElement]
    ): ErrorOr[Set[WomFile]] = {
      def filesInObjectField(fieldAndWomTypeTuple: (String, WomType)): ErrorOr[Set[WomFile]] = {
        val (field, womType) = fieldAndWomTypeTuple
        a.elements.get(field) match {
          case Some(fieldElement) =>
            fileEvaluator.predictFilesNeededToEvaluate(fieldElement, inputs, ioFunctionSet, womType)(fileEvaluator,
                                                                                                     valueEvaluator
            )
          case None => s"Invalid assignment to struct. Required field $field was not specified.".invalidNel
        }
      }

      coerceTo match {
        case WomCompositeType(mapping, _) => mapping.toList.traverse(filesInObjectField).map(_.flatten.toSet)
        case _ =>
          a.elements.values.toList
            .traverse(
              fileEvaluator.evaluateFilesNeededToEvaluate(_, inputs, ioFunctionSet, coerceTo)(fileEvaluator,
                                                                                              valueEvaluator
              )
            )
            .map(_.toSet.flatten)
      }
    }
  }
}
