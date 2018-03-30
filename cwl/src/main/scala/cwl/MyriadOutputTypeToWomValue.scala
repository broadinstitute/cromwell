package cwl

import cats.data.Validated.Valid
import cats.instances.list._
import cats.syntax.option._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import cwl.CwlType.CwlType
import cwl.MyriadOutputTypeToWomValue.EvaluationFunction
import cwl.command.ParentName
import mouse.all._
import shapeless.Poly1
import wom.types._
import wom.values.{WomObject, WomValue}

/**
  * Folds a MyriadOutputType into a WomValue
  * This is needed because the type might define the structure of the final output value (for OutputRecordSchemas for example)
  */
object MyriadOutputTypeToWomValue extends Poly1 {

  // We pass in a function that can evaluate a CommandOutputBinding and produce a WomValue. This allows us to recurse into the
  // MyriadOutputTypes and evaluate values as we do.
  type EvaluationFunction = (CommandOutputBinding, WomType) => ErrorOr[WomValue]

  import Case._

  implicit def cwlType: Aux[MyriadOutputInnerType, EvaluationFunction => ErrorOr[WomValue]] = at[MyriadOutputInnerType]{
    _.fold(MyriadOutputInnerTypeToWomValue)
  }

  // TODO: Not sure what the right thing to do is here, for now go over the list of types and use the first evaluation that yields success
  implicit def acwl: Aux[Array[MyriadOutputInnerType], EvaluationFunction => ErrorOr[WomValue]] = at[Array[MyriadOutputInnerType]] { types =>
    evalFunction =>
      types.toList.map(_.fold(MyriadOutputInnerTypeToWomValue).apply(evalFunction)).collectFirst({
        case Valid(validValue) => validValue
      }).toValidNel(s"Cannot find a suitable type to build a WomValue from in ${types.mkString(", ")}")
  }
}

object MyriadOutputInnerTypeToWomValue extends Poly1 {

  import Case._

  def ex(component: String) = throw new RuntimeException(s"output type $component cannot yield a wom value")

  implicit def cwlType: Aux[CwlType, EvaluationFunction => ErrorOr[WomValue]] = at[CwlType] { cwlType => _ =>
    "No output binding is defined. Are you expecting the output to be inferred from a cwl.output.json file ? If so please make sure the file was effectively created.".invalidNel
  }

  implicit def ors: Aux[OutputRecordSchema, EvaluationFunction => ErrorOr[WomValue]] = at[OutputRecordSchema] {
    case OutputRecordSchema(_, Some(fields), _) =>
      evalFunction =>

        // Go over each field and evaluate the binding if there's one, otherwise keep folding over field types
        def evaluateValues = fields.toList.traverse[ErrorOr, ((String, WomValue), (String, WomType))]({ field =>
          val womType = field.`type`.fold(MyriadOutputTypeToWomType)
          val womValue: ErrorOr[WomValue] = field.outputBinding match {
            case Some(binding) => evalFunction(binding, womType)
            case None => field.`type`.fold(MyriadOutputTypeToWomValue).apply(evalFunction)
          }

          // TODO: ParentName might need to be passed in here ?
          // return the value and the type with a clean parsedName
          womValue map { value =>
            val parsedName = FullyQualifiedName(field.name)(ParentName.empty).id
            (parsedName -> value) -> (parsedName -> womType)
          }
        })

        evaluateValues map { evaluatedValues =>
          val (valueMap, typeMap) = evaluatedValues.unzip
          // Create a typed WomObject from the values and the typeMap
          WomObject.withType(valueMap.toMap, WomCompositeType(typeMap.toMap))
        }
    case ors => ors.toString |> ex
  }

  implicit def oes: Aux[OutputEnumSchema, EvaluationFunction => ErrorOr[WomValue]] = at[OutputEnumSchema]{
    //DB: I tried to do a pattern match as the overall function here but the compiler exploded
    oes => f => oes match {
      case oes@OutputEnumSchema(_, _, _, Some(outputBinding)) => f(outputBinding, oes.toWomEnumerationType)
      case _ => s"The enumeration type $oes requires an outputbinding to be evaluated.".invalidNel
    }
  }

  implicit def oas: Aux[OutputArraySchema, EvaluationFunction => ErrorOr[WomValue]] = at[OutputArraySchema]{
    case OutputArraySchema(itemsType, _, _, outputBinding) =>
      evalFunction =>
        lazy val itemsWomType = itemsType.fold(MyriadOutputTypeToWomType)
        def fromBinding = outputBinding.map(evalFunction(_, WomArrayType(itemsWomType)))
        def fromTypes = itemsType.fold(MyriadOutputTypeToWomValue).apply(evalFunction)
        fromBinding.getOrElse(fromTypes)
  }

  implicit def s: Aux[String, EvaluationFunction => ErrorOr[WomValue]] = at[String]{ s => _ =>
    s.toString |> ex
  }
}
