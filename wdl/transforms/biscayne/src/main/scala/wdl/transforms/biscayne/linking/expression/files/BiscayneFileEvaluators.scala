package wdl.transforms.biscayne.linking.expression.files

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
import wdl.transforms.base.linking.expression.files.EngineFunctionEvaluators.{
  threeParameterFunctionPassthroughFileEvaluator,
  twoParameterFunctionPassthroughFileEvaluator
}
import wdl.transforms.base.linking.expression.files.EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
import wom.expression.IoFunctionSet
import wom.types.{WomCompositeType, WomType}
import wom.values.{WomFile, WomValue}

object BiscayneFileEvaluators {

  implicit val keysFileEvaluator: FileEvaluator[Keys] = singleParameterPassthroughFileEvaluator
  implicit val asMapFileEvaluator: FileEvaluator[AsMap] = singleParameterPassthroughFileEvaluator
  implicit val asPairsFileEvaluator: FileEvaluator[AsPairs] = singleParameterPassthroughFileEvaluator
  implicit val collectByKeyFileEvaluator: FileEvaluator[CollectByKey] = singleParameterPassthroughFileEvaluator

  implicit val sepFunctionEvaluator: FileEvaluator[Sep] = twoParameterFunctionPassthroughFileEvaluator[Sep]
  implicit val subPosixFunctionEvaluator: FileEvaluator[SubPosix] =
    threeParameterFunctionPassthroughFileEvaluator[SubPosix]
  implicit val suffixFunctionEvaluator: FileEvaluator[Suffix] = twoParameterFunctionPassthroughFileEvaluator[Suffix]
  implicit val quoteFunctionEvaluator: FileEvaluator[Quote] = singleParameterPassthroughFileEvaluator
  implicit val sQuoteFunctionEvaluator: FileEvaluator[SQuote] = singleParameterPassthroughFileEvaluator

  implicit val minFunctionEvaluator: FileEvaluator[Min] = twoParameterFunctionPassthroughFileEvaluator[Min]
  implicit val maxFunctionEvaluator: FileEvaluator[Max] = twoParameterFunctionPassthroughFileEvaluator[Max]

  implicit val unzipFunctionEvaluator: FileEvaluator[Unzip] = singleParameterPassthroughFileEvaluator

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
