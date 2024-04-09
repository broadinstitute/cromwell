package wdl.transforms.biscayne.linking.expression.types

import cats.data.Validated.{Invalid, Valid}
import cats.implicits.{catsSyntaxTuple2Semigroupal, catsSyntaxTuple3Semigroupal}
import cats.syntax.validated._
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.transforms.base.linking.expression.types.EngineFunctionEvaluators.validateParamType
import wom.types._

object BiscayneTypeEvaluators {
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
          s"Cannot invoke squote on type Array[${other.stableName}]. Expected an Array of primitive type".invalidNel
        case other =>
          s"Cannot invoke squote on type ${other.stableName}. Expected an Array of primitive type".invalidNel
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

    // does it make sense that someone would assign type b to type a?
    def areTypesAssignable(a: WomType, b: WomType): Boolean =
      !a.equalsType(b).isFailure

    // Helper method to check something (maybe) found in the struct literal to something (maybe) found in the struct definition.
    def checkIfMemberIsValid(typeName: String,
                             memberName: String,
                             evaluatedType: Option[WomType],
                             expectedType: Option[WomType]
    ): ErrorOr[WomType] =
      evaluatedType match {
        case Some(evaluated) =>
          expectedType match {
            case Some(expected) =>
              if (areTypesAssignable(evaluated, expected)) evaluated.validNel
              else
                s"$typeName.$memberName expected to be ${expected.friendlyName}. Found ${evaluated.friendlyName}.".invalidNel
            case None => s"Type $typeName does not have a member called $memberName.".invalidNel
          }
        case None => s"Error evaluating the type of ${typeName}.${memberName}.".invalidNel
      }

    // For each member in the literal, check that it exists in the struct definition and is the expected type.
    def checkMembersAgainstDefinition(a: StructLiteral,
                                      structDefinition: WomCompositeType,
                                      linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                      typeAliases: Map[String, WomType]
    ): ErrorOr[WomCompositeType] = {
      val checkedMembers: Map[String, ErrorOr[WomType]] = a.elements.map { case (memberKey, memberExpressionElement) =>
        val evaluatedType =
          expressionTypeEvaluator.evaluateType(memberExpressionElement, linkedValues, typeAliases).toOption
        val expectedType = structDefinition.typeMap.get(memberKey)
        (memberKey, checkIfMemberIsValid(a.structTypeName, memberKey, evaluatedType, expectedType))
      }

      val errors: Iterable[String] = checkedMembers.flatMap { case (_, errorOr) =>
        errorOr match {
          case Invalid(e) => Some(e.toList.mkString)
          case _ => None
        }
      }

      if (errors.nonEmpty) {
        errors.mkString(",").invalidNel
      } else {
        val validatedTypes: Map[String, WomType] = checkedMembers.flatMap { case (key, errorOr) =>
          errorOr match {
            case Valid(v) => Some((key, v))
            case _ => None
          }
        }
        WomCompositeType(validatedTypes, Some(a.structTypeName)).validNel
      }
    }

    // For every member in the definition, if that member isn't optional, confirm that it is also in the struct literal.
    def checkForMissingMembers(foundMembers: Map[String, WomType],
                               structDefinition: WomCompositeType
    ): Option[String] = {
      val errors: Iterable[String] = structDefinition.typeMap flatMap { case (memberName, memberType) =>
        memberType match {
          case WomOptionalType(_) => None
          case _ =>
            if (!foundMembers.contains(memberName)) Some(s"Expected member ${memberName} not found. ")
            else None
        }
      }
      errors match {
        case Nil => None
        case _ => Some(errors.mkString)
      }
    }
    override def evaluateType(a: StructLiteral,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] = {
      val structDefinition = typeAliases.get(a.structTypeName)
      structDefinition match {
        case Some(definition) =>
          definition match {
            case compositeType: WomCompositeType =>
              checkMembersAgainstDefinition(a, compositeType, linkedValues, typeAliases).flatMap { foundMembers =>
                checkForMissingMembers(foundMembers.typeMap, compositeType) match {
                  case Some(error) => error.invalidNel
                  case _ => compositeType.validNel
                }
              }
            case _ => s"Unexpected error while parsing ${a.structTypeName}".invalidNel
          }
        case None => s"Could not find Struct Definition for type ${a.structTypeName}".invalidNel
      }
    }
  }
}
