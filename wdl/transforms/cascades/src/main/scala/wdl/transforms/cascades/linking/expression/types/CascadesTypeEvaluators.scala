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

object cascadesTypeEvaluators {
  implicit val keysFunctionEvaluator: TypeEvaluator[Keys] = new TypeEvaluator[Keys] {
    override def evaluateType(a: Keys,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      validateParamType(a.param, linkedValues, WomMapType(WomAnyType, WomAnyType), typeAliases) flatMap {
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
      validateParamType(a.param, linkedValues, WomArrayType(WomPairType(WomAnyType, WomAnyType)), typeAliases) flatMap {
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
      validateParamType(a.param, linkedValues, WomMapType(WomAnyType, WomAnyType), typeAliases) flatMap {
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
      validateParamType(a.param, linkedValues, WomArrayType(WomPairType(WomAnyType, WomAnyType)), typeAliases) flatMap {
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
      validateParamType(a.arg2, linkedValues, WomArrayType(WomAnyType), typeAliases) flatMap {
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
      (validateParamType(a.input, linkedValues, WomSingleFileType, typeAliases),
       validateParamType(a.pattern, linkedValues, WomSingleFileType, typeAliases),
       validateParamType(a.replace, linkedValues, WomSingleFileType, typeAliases)
      ) mapN { (_, _, _) => WomStringType }
  }

  implicit val suffixFunctionEvaluator: TypeEvaluator[Suffix] = new TypeEvaluator[Suffix] {
    override def evaluateType(a: Suffix,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      (validateParamType(a.suffix, linkedValues, WomStringType, typeAliases),
       validateParamType(a.array, linkedValues, WomArrayType(WomStringType), typeAliases)
      ) mapN { (_, _) => WomArrayType(WomStringType) }
  }

  implicit val quoteFunctionEvaluator: TypeEvaluator[Quote] = new TypeEvaluator[Quote] {
    override def evaluateType(a: Quote,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      validateParamType(a.param, linkedValues, WomArrayType(WomAnyType), typeAliases) flatMap {
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
      validateParamType(a.param, linkedValues, WomArrayType(WomAnyType), typeAliases) flatMap {
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
      validateParamType(a.param, linkedValues, WomArrayType(WomPairType(WomAnyType, WomAnyType)), typeAliases) flatMap {
        case WomArrayType(WomNothingType) => WomPairType(WomArrayType(WomAnyType), WomArrayType(WomAnyType)).validNel
        case WomArrayType(WomPairType(x, y)) => WomPairType(WomArrayType(x), WomArrayType(y)).validNel
        case other => s"Cannot invoke 'unzip' on type '${other.stableName}'. Expected an array of pairs".invalidNel
      }
  }

  implicit val structLiteralTypeEvaluator: TypeEvaluator[StructLiteral] = new TypeEvaluator[StructLiteral] {

    /**
     * Is the evaluated type allowed to be assigned to the expectedType?
     */
    def areTypesAssignable(evaluatedType: WomType, expectedType: WomType): Boolean =
      // NB: This check is a little looser than we'd like it to be.
      // For example, String is coercible to Int (Int i = "1" is OK)
      // It's not until we actually evaluate the value of the string that we can know if that coercion succeeded or not. (Int i = "orange" will fail)
      // We don't know whether the user has provided "1" or "orange" at this stage.
      // This is OK as-is because the value evaluators do the coercing and throw meaningful errors if the coercion fails.
      expectedType.isCoerceableFrom(evaluatedType)

    /**
     * Helper method to check if something (maybe) found in the struct literal against something (maybe) found in the struct definition.
     * @return The WomType of the evaluated member. Error if either is not present, or if the types aren't compatible.
     */
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

    /**
     * For each member in the literal, check that it exists in the struct definition and is the expected type.
     * @return The WomCompositeType of the struct literal, as determined from evaluating each member.
     *         This might not *exactly* match the struct definition due to permitted type coercions.
     */
    def checkMembersAgainstDefinition(a: StructLiteral,
                                      structDefinition: WomCompositeType,
                                      linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                      typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomCompositeType] = {
      val checkedMembers: Map[String, ErrorOr[WomType]] = a.elements.map { case (memberKey, memberExpressionElement) =>
        val evaluatedType =
          expressionTypeEvaluator.evaluateType(memberExpressionElement, linkedValues, typeAliases).toOption
        val expectedType = structDefinition.typeMap.get(memberKey)
        (memberKey, checkIfMemberIsValid(a.structTypeName, memberKey, evaluatedType, expectedType))
      }

      val (errors, validatedTypes) = checkedMembers.partition { case (_, errorOr) =>
        errorOr match {
          case Invalid(_) => true
          case Valid(_) => false
        }
      }

      if (errors.nonEmpty) {
        errors.collect { case (_, Invalid(e)) => e.toList.mkString(", ") }.toList.mkString("[ ", ", ", " ]").invalidNel
      } else {
        val types = validatedTypes.collect { case (key, Valid(v)) => (key, v) }
        WomCompositeType(types, Some(a.structTypeName)).validNel
      }
    }

    /**
     * For each member in the struct definition, if that member isn't optional, confirm that it is also in the struct literal.
     */
    def checkForMissingMembers(foundMembers: Map[String, WomType],
                               structDefinition: WomCompositeType
    ): ErrorOr[Unit] = {
      val errors: Iterable[String] = structDefinition.typeMap flatMap { case (memberName, memberType) =>
        memberType match {
          case WomOptionalType(_) => None
          case _ =>
            if (!foundMembers.contains(memberName)) Some(s"Expected member ${memberName} not found. ")
            else None
        }
      }
      if (errors.nonEmpty) errors.mkString.invalidNel else ().validNel
    }

    /**
     * Returns the type of a struct literal, assuming it is compatible with an existing struct definition.
     * It is required that:
     *   - a struct definition exist for the literal (looked up the via name of the struct type). This is defined elsewhere in the WDL and is not part of the literal.
     *   - the struct definition be a WomCompositeType (it's programmer error if it's not)
     *   - each  member provided in the literal is also in the definition, and is coercible to the defined type.
     *   - all non-optional members of the definition are present in the literal.
     * @param a The literal to evaluate
     * @param linkedValues Used by the expression type evaluator.
     * @param typeAliases A map containing available struct definitions
     * @param expressionTypeEvaluator An object capable of evaluating the types of ExpressionElements. Used to evaluate the type of each member provided in the literal.
     * @return The type of the struct definition (as found in typeAliases) if all goes well. A list of errors otherwise.
     */
    override def evaluateType(a: StructLiteral,
                              linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                              typeAliases: Map[String, WomType]
    )(implicit
      expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
    ): ErrorOr[WomType] =
      typeAliases.get(a.structTypeName) match {
        case Some(definition) =>
          definition match {
            case compositeType: WomCompositeType =>
              checkMembersAgainstDefinition(a, compositeType, linkedValues, typeAliases).flatMap { foundMembers =>
                checkForMissingMembers(foundMembers.typeMap, compositeType) match {
                  case Invalid(error) => error.invalid
                  case _ => compositeType.validNel
                }
              }
            case _ =>
              s"Programmer error: Expected the struct definition of ${a.structTypeName} to be a WomCompositeType".invalidNel
          }
        case None => s"Could not find Struct Definition for type ${a.structTypeName}".invalidNel
      }
  }
}
