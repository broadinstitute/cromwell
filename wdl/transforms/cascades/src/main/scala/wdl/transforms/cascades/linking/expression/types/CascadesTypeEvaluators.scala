package wdl.transforms.biscayne.linking.expression.types

import cats.implicits.{catsSyntaxTuple2Semigroupal, catsSyntaxTuple3Semigroupal}
import cats.syntax.validated._
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.transforms.base.linking.expression.types.EngineFunctionEvaluators.validateParamType
import wom.types._

object cascadesTypeEvaluators {
  implicit val keysFunctionEvaluator: TypeEvaluator[Keys] = new TypeEvaluator[Keys] {
    override def evaluateType(a: Keys,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      validateParamType(a.param, linkedValues, WomMapType(WomAnyType, WomAnyType)) flatMap {
        case WomMapType(keyType, _) => WomArrayType(keyType).validNel
        case other => s"Cannot invoke 'keys' on type '${other.stableName}'. Expected a map".invalidNel
      }
  }

  implicit val asMapFunctionEvaluator: TypeEvaluator[AsMap] = new TypeEvaluator[AsMap] {
    override def evaluateType(a: AsMap,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      validateParamType(a.param, linkedValues, WomArrayType(WomPairType(WomAnyType, WomAnyType))) flatMap {
        case WomArrayType(WomPairType(x: WomPrimitiveType, y)) => WomMapType(x, y).validNel
        case other @ WomArrayType(WomPairType(x, _)) =>
          s"Cannot invoke 'as_map' on type ${other.stableName}. Map keys must be primitive but got '${x.stableName}'".invalidNel
        case other => s"Cannot invoke 'as_map' on type '${other.stableName}'. Expected an array of pairs".invalidNel
      }
  }

  implicit val asPairsFunctionEvaluator: TypeEvaluator[AsPairs] = new TypeEvaluator[AsPairs] {
    override def evaluateType(a: AsPairs,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      validateParamType(a.param, linkedValues, WomMapType(WomAnyType, WomAnyType)) flatMap {
        case WomMapType(x, y) => WomArrayType(WomPairType(x, y)).validNel
        case other => s"Cannot invoke 'as_pairs' on type '${other.stableName}'. Expected a map".invalidNel
      }
  }

  implicit val collectByKeyFunctionEvaluator: TypeEvaluator[CollectByKey] = new TypeEvaluator[CollectByKey] {
    override def evaluateType(a: CollectByKey,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      validateParamType(a.param, linkedValues, WomArrayType(WomPairType(WomAnyType, WomAnyType))) flatMap {
        case WomArrayType(WomPairType(x: WomPrimitiveType, y)) => WomMapType(x, WomArrayType(y)).validNel
        case other @ WomArrayType(WomPairType(x, _)) =>
          s"Cannot invoke 'collect_by_key' on type ${other.stableName}. Map keys must be primitive but got '${x.stableName}'".invalidNel
        case other =>
          s"Cannot invoke 'collect_by_key' on type '${other.stableName}'. Expected an array of pairs".invalidNel
      }
  }

  private def resultTypeOfIntVsFloat(functionName: String)(type1: WomType, type2: WomType): ErrorOr[WomType] =
    (type1, type2) match {
      case (WomIntegerType, WomIntegerType) => WomIntegerType.validNel
      case (WomIntegerType, WomFloatType) => WomFloatType.validNel
      case (WomFloatType, WomIntegerType) => WomFloatType.validNel
      case (WomFloatType, WomFloatType) => WomFloatType.validNel
      case (other1, other2) =>
        s"Cannot call '$functionName' with arguments (${other1.friendlyName}, ${other2.friendlyName}). Must be Int or Long.".invalidNel
    }

  implicit val minFunctionEvaluator: TypeEvaluator[Min] = new TypeEvaluator[Min] {
    override def evaluateType(a: Min,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] = {
      val type1 = expressionTypeEvaluator.evaluateType(a.arg1, linkedValues, typeAliases)
      val type2 = expressionTypeEvaluator.evaluateType(a.arg1, linkedValues, typeAliases)

      (type1, type2) flatMapN resultTypeOfIntVsFloat("min")
    }
  }

  implicit val maxFunctionEvaluator: TypeEvaluator[Max] = new TypeEvaluator[Max] {
    override def evaluateType(a: Max,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] = {
      val type1 = expressionTypeEvaluator.evaluateType(a.arg1, linkedValues, typeAliases)
      val type2 = expressionTypeEvaluator.evaluateType(a.arg1, linkedValues, typeAliases)

      (type1, type2) flatMapN resultTypeOfIntVsFloat("max")
    }
  }

  implicit val sepFunctionEvaluator: TypeEvaluator[Sep] = new TypeEvaluator[Sep] {
    override def evaluateType(a: Sep,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      validateParamType(a.arg2, linkedValues, WomArrayType(WomAnyType)) flatMap {
        case WomArrayType(WomArrayType(_)) =>
          s"Cannot invoke 'sep' on type 'Array[Array[_]]'. Expected an Array[String].".invalidNel
        case WomArrayType(_) => WomStringType.validNel

        case other => s"Cannot invoke 'sep' on type '${other.stableName}'. Expected an Array[String].".invalidNel

      }
  }

  implicit val subPosixFunctionEvaluator: TypeEvaluator[SubPosix] = new TypeEvaluator[SubPosix] {
    override def evaluateType(a: SubPosix,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      (validateParamType(a.input, linkedValues, WomSingleFileType),
       validateParamType(a.pattern, linkedValues, WomSingleFileType),
       validateParamType(a.replace, linkedValues, WomSingleFileType)
      ) mapN { (_, _, _) => WomStringType }
  }

  implicit val suffixFunctionEvaluator: TypeEvaluator[Suffix] = new TypeEvaluator[Suffix] {
    override def evaluateType(a: Suffix,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      (validateParamType(a.suffix, linkedValues, WomStringType),
       validateParamType(a.array, linkedValues, WomArrayType(WomStringType))
      ) mapN { (_, _) => WomArrayType(WomStringType) }
  }

  implicit val quoteFunctionEvaluator: TypeEvaluator[Quote] = new TypeEvaluator[Quote] {
    override def evaluateType(a: Quote,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      validateParamType(a.param, linkedValues, WomArrayType(WomAnyType)) flatMap {
        case WomArrayType(WomNothingType) => WomArrayType(WomNothingType).validNel
        case WomArrayType(_: WomPrimitiveType) => WomArrayType(WomStringType).validNel
        case other @ WomArrayType(_) =>
          s"Cannot invoke quote on type Array[${other.stableName}]. Expected an Array of primitive type".invalidNel
        case other =>
          s"Cannot invoke quote on type ${other.stableName}. Expected an Array of primitive type".invalidNel
      }
  }

  implicit val sQuoteFunctionEvaluator: TypeEvaluator[SQuote] = new TypeEvaluator[SQuote] {
    override def evaluateType(a: SQuote,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      validateParamType(a.param, linkedValues, WomArrayType(WomAnyType)) flatMap {
        case WomArrayType(WomNothingType) => WomArrayType(WomNothingType).validNel
        case WomArrayType(_: WomPrimitiveType) => WomArrayType(WomStringType).validNel
        case other @ WomArrayType(_) =>
          s"Cannot invoke quote on type Array[${other.stableName}]. Expected an Array of primitive type".invalidNel
        case other =>
          s"Cannot invoke quote on type ${other.stableName}. Expected an Array of primitive type".invalidNel
      }
  }

  implicit val unzipFunctionEvaluator: TypeEvaluator[Unzip] = new TypeEvaluator[Unzip] {
    override def evaluateType(a: Unzip,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      validateParamType(a.param, linkedValues, WomArrayType(WomPairType(WomAnyType, WomAnyType))) flatMap {
        case WomArrayType(WomNothingType) => WomPairType(WomArrayType(WomAnyType), WomArrayType(WomAnyType)).validNel
        case WomArrayType(WomPairType(x, y)) => WomPairType(WomArrayType(x), WomArrayType(y)).validNel
        case other => s"Cannot invoke 'unzip' on type '${other.stableName}'. Expected an array of pairs".invalidNel
      }
  }

  implicit val structLiteralTypeEvaluator: TypeEvaluator[StructLiteral] = new TypeEvaluator[StructLiteral] {
    override def evaluateType(a: StructLiteral,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      // This works fine, but is not yet a strong enough type check for the WDL 1.1 spec
      // (i.e. users are able to instantiate struct literals with k/v pairs that aren't actually in the struct definition, without an error being thrown.)
      // We want to add extra validation here, and return a WomCompositeType that matches the struct definition of everything is OK.
      // Note that users are allowed to omit optional k/v pairs in their literal.
      // This requires some extra work to be done in a subsequent PR.
      WomObjectType.validNel
  }
}
