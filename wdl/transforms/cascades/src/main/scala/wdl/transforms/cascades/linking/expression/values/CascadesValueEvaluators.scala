package wdl.transforms.biscayne.linking.expression.values

import cats.data.NonEmptyList
import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import com.google.re2j.{Pattern => RE2JPattern}
import common.validation.ErrorOr._
import common.collections.EnhancedCollections._
import common.validation.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{EvaluatedValue, ForCommandInstantiationOptions, ValueEvaluator}
import wdl.transforms.base.linking.expression.values.EngineFunctionEvaluators.{
  processThreeValidatedValues,
  processTwoValidatedValues,
  processValidatedSingleValue
}
import wom.expression.IoFunctionSet
import wom.types._
import wom.values.{WomArray, WomFloat, WomInteger, WomMap, WomObject, WomOptionalValue, WomPair, WomString, WomValue}
import wom.types.coercion.defaults._

object cascadesValueEvaluators {

  implicit val noneLiteralEvaluator: ValueEvaluator[NoneLiteralElement.type] =
    new ValueEvaluator[ExpressionElement.NoneLiteralElement.type] {
      override def evaluateValue(a: ExpressionElement.NoneLiteralElement.type,
                                 inputs: Map[String, WomValue],
                                 ioFunctionSet: IoFunctionSet,
                                 forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
      )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] =
        EvaluatedValue(value = WomOptionalValue(WomNothingType, None), sideEffectFiles = Seq.empty).validNel
    }

  implicit val asMapFunctionEvaluator: ValueEvaluator[AsMap] = new ValueEvaluator[AsMap] {
    override def evaluateValue(a: AsMap,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] =
      processValidatedSingleValue[WomArray, WomMap](
        expressionValueEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )
      ) {
        case WomArray(WomArrayType(WomPairType(_: WomPrimitiveType, _)), values) =>
          val validPairs: ErrorOr[List[(WomValue, WomValue)]] = values.toList traverse {
            case WomPair(l, r) => (l, r).validNel
            case other =>
              s"Unexpected array element. Expected a Pair[X, Y] but array contained ${other.toWomString}]".invalidNel
          } leftMap { errors =>
            NonEmptyList.fromListUnsafe(errors.toList.distinct)
          }
          validPairs flatMap { pairs =>
            val grouped = pairs.groupBy(_._1)
            val tooManyKeyErrors = grouped collect {
              case (name, list) if list.length != 1 =>
                s"keys can only appear once but ${name.toWomString} appeared ${list.size} times."
            }
            if (tooManyKeyErrors.isEmpty) {
              val pairs = grouped map { case (key, value) => key -> value.head._2 }

              EvaluatedValue(WomMap(pairs), Seq.empty).validNel
            } else {
              s"Cannot evaluate 'as_map' with duplicated keys: ${tooManyKeyErrors.mkString(", ")}".invalidNel
            }
          }

        case WomArray(womType @ WomArrayType(WomPairType(x, _)), _) =>
          s"Cannot evaluate 'as_map' on type ${womType.stableName}. Keys must be primitive but got ${x.stableName}.".invalidNel
        case other =>
          s"Invalid call of 'as_map' on parameter of type '${other.womType.stableName}' (expected Array[Pair[X, Y]])".invalidNel
      }(coercer = WomArrayType(WomPairType(WomAnyType, WomAnyType)))
  }

  implicit val keysFunctionEvaluator: ValueEvaluator[Keys] = new ValueEvaluator[Keys] {
    override def evaluateValue(a: Keys,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] =
      processValidatedSingleValue[WomMap, WomArray](
        expressionValueEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )
      ) {
        case WomMap(WomMapType(keyType, _), values) =>
          EvaluatedValue(WomArray(WomArrayType(keyType), values.keys.toList), Seq.empty).validNel
        case other =>
          s"Invalid call of 'keys' on parameter of type '${other.womType.stableName}' (expected Map[X, Y])".invalidNel
      }
  }

  implicit val asPairsFunctionEvaluator: ValueEvaluator[AsPairs] = new ValueEvaluator[AsPairs] {
    override def evaluateValue(a: AsPairs,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] =
      processValidatedSingleValue[WomMap, WomArray](
        expressionValueEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )
      ) {
        case WomMap(WomMapType(keyType, valueType), values) =>
          val validPairs: List[WomPair] = values.toList map { case (l, r) =>
            WomPair(l, r)
          }
          EvaluatedValue(WomArray(WomArrayType(WomPairType(keyType, valueType)), validPairs), Seq.empty).validNel

        case other =>
          s"Invalid call of 'as_pairs' on parameter of type '${other.womType.stableName}' (expected Map[X, Y])".invalidNel
      }
  }

  implicit val collectByKeyFunctionEvaluator: ValueEvaluator[CollectByKey] = new ValueEvaluator[CollectByKey] {
    override def evaluateValue(a: CollectByKey,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] =
      processValidatedSingleValue[WomArray, WomMap](
        expressionValueEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )
      ) {
        case WomArray(WomArrayType(WomPairType(_: WomPrimitiveType, _)), values) =>
          val validPairs: ErrorOr[List[(WomValue, WomValue)]] = values.toList traverse {
            case WomPair(l, r) => (l, r).validNel
            case other =>
              s"Unexpected array element. Expected a Pair[X, Y] but array contained ${other.toWomString}]".invalidNel
          }
          validPairs flatMap { kvpairs =>
            val grouped: Map[WomValue, WomArray] = kvpairs.groupBy(_._1).safeMapValues(v => WomArray(v.map(_._2)))
            EvaluatedValue(WomMap(grouped), Seq.empty).validNel

          }

        case WomArray(womType @ WomArrayType(WomPairType(x, _)), _) =>
          s"Cannot evaluate 'collect_by_key' on type ${womType.stableName}. Keys must be primitive but got ${x.stableName}.".invalidNel
        case other =>
          s"Invalid call of 'collect_by_key' on parameter of type '${other.womType.stableName}' (expected Map[X, Y])".invalidNel
      }(coercer = WomArrayType(WomPairType(WomAnyType, WomAnyType)))
  }

  private def resultOfIntVsFloat(functionName: String,
                                 intFunc: (Int, Int) => Int,
                                 doubleFunc: (Double, Double) => Double
  )(value1: EvaluatedValue[_], value2: EvaluatedValue[_]): ErrorOr[EvaluatedValue[WomValue]] = {
    val newValue = (value1.value, value2.value) match {
      case (WomInteger(i1), WomInteger(i2)) => WomInteger(intFunc(i1, i2)).validNel
      case (WomInteger(i1), WomFloat(l2)) => WomFloat(doubleFunc(i1.doubleValue, l2)).validNel
      case (WomFloat(l1), WomInteger(i2)) => WomFloat(doubleFunc(l1, i2.doubleValue)).validNel
      case (WomFloat(l1), WomFloat(l2)) => WomFloat(doubleFunc(l1, l2)).validNel
      case (other1, other2) =>
        s"Invalid arguments to '$functionName':(${other1.typeName}, ${other2.typeName})".invalidNel
    }
    newValue map { v => EvaluatedValue(v, value1.sideEffectFiles ++ value2.sideEffectFiles) }
  }

  implicit val minFunctionEvaluator: ValueEvaluator[Min] = new ValueEvaluator[Min] {
    override def evaluateValue(a: Min,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomValue]] = {
      val value1 =
        expressionValueEvaluator.evaluateValue(a.arg1, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )
      val value2 =
        expressionValueEvaluator.evaluateValue(a.arg2, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )

      val intFunc = (i1: Int, i2: Int) => Math.min(i1, i2)
      val doubleFunc = (l1: Double, l2: Double) => Math.min(l1, l2)

      (value1, value2) flatMapN resultOfIntVsFloat("min", intFunc, doubleFunc)
    }
  }

  implicit val maxFunctionEvaluator: ValueEvaluator[Max] = new ValueEvaluator[Max] {
    override def evaluateValue(a: Max,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomValue]] = {
      val value1 =
        expressionValueEvaluator.evaluateValue(a.arg1, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )
      val value2 =
        expressionValueEvaluator.evaluateValue(a.arg2, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )

      val intFunc = (i1: Int, i2: Int) => Math.max(i1, i2)
      val doubleFunc = (l1: Double, l2: Double) => Math.max(l1, l2)

      (value1, value2) flatMapN resultOfIntVsFloat("max", intFunc, doubleFunc)
    }
  }

  implicit val sepFunctionEvaluator: ValueEvaluator[Sep] = new ValueEvaluator[Sep] {
    override def evaluateValue(a: Sep,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomString]] =
      processTwoValidatedValues[WomString, WomArray, WomString](
        expressionValueEvaluator.evaluateValue(a.arg1, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        ),
        expressionValueEvaluator.evaluateValue(a.arg2, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )
      ) { (sepvalue, arr) =>
        EvaluatedValue(WomString(arr.value.map(v => v.valueString).mkString(sepvalue.value)), Seq.empty).validNel
      }
  }

  implicit val subPosixFunctionEvaluator: ValueEvaluator[SubPosix] = new ValueEvaluator[SubPosix] {
    override def evaluateValue(a: SubPosix,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomString]] =
      processThreeValidatedValues[WomString, WomString, WomString, WomString](
        expressionValueEvaluator.evaluateValue(a.input, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        ),
        expressionValueEvaluator.evaluateValue(a.pattern, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        ),
        expressionValueEvaluator.evaluateValue(a.replace, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )
      ) { (input, pattern, replace) =>
        ErrorOr(
          EvaluatedValue(WomString(
                           RE2JPattern
                             .compile(pattern.valueString)
                             .matcher(input.valueString)
                             .replaceAll(replace.valueString)
                         ),
                         Seq.empty
          )
        )
      }
  }

  implicit val suffixFunctionEvaluator: ValueEvaluator[Suffix] = new ValueEvaluator[Suffix] {
    override def evaluateValue(a: Suffix,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] =
      processTwoValidatedValues[WomString, WomArray, WomArray](
        expressionValueEvaluator.evaluateValue(a.arg1, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        ),
        expressionValueEvaluator.evaluateValue(a.arg2, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )
      ) { (suffix, arr) =>
        EvaluatedValue(WomArray(arr.value.map(v => WomString(v.valueString + suffix.value))), Seq.empty).validNel
      }
  }

  /**
   * Quote: Given an array of primitives, produce a new array in which all elements of the original are in quotes (").
   * https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-arraystring-quotearrayp
   * input: Array[Primitive]
   * output: Array[String]
   */
  implicit val quoteFunctionEvaluator: ValueEvaluator[Quote] = new ValueEvaluator[Quote] {
    override def evaluateValue(a: Quote,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] =
      processValidatedSingleValue[WomArray, WomArray](
        expressionValueEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )
      ) { arr =>
        EvaluatedValue(WomArray(arr.value.map(v => WomString("\"" + v.valueString + "\""))), Seq.empty).validNel
      }
  }

  /**
   * SQuote: Given an array of primitives, produce a new array in which all elements of the original are in single quotes (').
   * https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-arraystring-squotearrayp
   * input: Array[Primitive]
   * output: Array[String]
   */
  implicit val sQuoteFunctionEvaluator: ValueEvaluator[SQuote] = new ValueEvaluator[SQuote] {
    override def evaluateValue(a: SQuote,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] =
      processValidatedSingleValue[WomArray, WomArray](
        expressionValueEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )
      ) { arr =>
        EvaluatedValue(WomArray(arr.value.map(v => WomString("\'" + v.valueString + "\'"))), Seq.empty).validNel
      }
  }

  /**
   * Unzip: Creates a pair of arrays, the first containing the elements from the left members of an array of pairs,
   * and the second containing the right members. This is the inverse of the zip function.
   * https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#-pairarrayx-arrayy-unziparraypairx-y
   * input: Array[Pair[X,Y]]
   * output: Pair[Array[X], Array[Y]]
   */
  implicit val unzipFunctionEvaluator: ValueEvaluator[Unzip] = new ValueEvaluator[Unzip] {
    override def evaluateValue(a: Unzip,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomPair]] =
      processValidatedSingleValue[WomArray, WomPair](
        expressionValueEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(
          expressionValueEvaluator
        )
      ) {
        case WomArray(WomArrayType(WomAnyType), Seq()) =>
          EvaluatedValue(
            WomPair(WomArray(WomArrayType(WomAnyType), Seq.empty), WomArray(WomArrayType(WomAnyType), Seq.empty)),
            Seq.empty
          ).validNel
        case WomArray(WomArrayType(WomPairType(_, _)), values) =>
          val zippedArrayOfPairs: Seq[(WomValue, WomValue)] = values map { case pair: WomPair =>
            Tuple2(pair.left, pair.right)
          }
          val (left, right) = zippedArrayOfPairs.unzip
          val unzippedPairOfArrays: WomPair = WomPair(WomArray(left), WomArray(right))
          EvaluatedValue(unzippedPairOfArrays, Seq.empty).validNel
        case other =>
          s"Invalid call of 'unzip' on parameter of type '${other.womType.stableName}' (expected Array[Pair[X, Y]])".invalidNel
      }
  }

  implicit val structLiteralValueEvaluator: ValueEvaluator[StructLiteral] = new ValueEvaluator[StructLiteral] {
    // This works fine, but is missing a feature from the WDL 1.1 spec: users are allowed to omit optional values from their struct literal.
    // This requires some extra work to be done in a subsequent PR.
    // Specifically, make the known struct definitions available to this function so we can populate k/v pairs with None appropriately.
    override def evaluateValue(a: StructLiteral,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
    )(implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {

      val evaluated: ErrorOr[List[(String, EvaluatedValue[_])]] = a.elements.toList traverse {
        case (key: String, value: ExpressionElement) =>
          expressionValueEvaluator
            .evaluateValue(value, inputs, ioFunctionSet, forCommandInstantiationOptions)(expressionValueEvaluator)
            .map(key -> _)
      }

      evaluated map { mapping =>
        val value = mapping.map(entry => entry._1 -> entry._2.value).toMap
        val sideEffectFiles = mapping.flatMap(entry => entry._2.sideEffectFiles)
        EvaluatedValue(WomObject(value), sideEffectFiles)
      }
    }
  }
}
