package cwl

import cats.data.Validated.Valid
import cats.instances.list._
import cats.syntax.option._
import cats.syntax.traverse._
import common.validation.IOChecked.IOChecked
import common.validation.IOChecked._
import cwl.CwlType.CwlType
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
  type EvaluationFunction = (CommandOutputBinding, WomType) => IOChecked[WomValue]

  //Our overall return type gives us the evaluator function and custom types; returns WomValues
  type Output = (EvaluationFunction, SchemaDefRequirement) => IOChecked[WomValue]

  import Case._

  implicit def cwlType: Aux[MyriadOutputInnerType, Output] = at[MyriadOutputInnerType]{
    moit => (evalFunction, schemaDefRequirement) => moit.fold(MyriadOutputInnerTypeToWomValue).apply(evalFunction, schemaDefRequirement)
  }

  // TODO: Not sure what the right thing to do is here, for now go over the list of types and use the first evaluation that yields success
  implicit def acwl: Aux[Array[MyriadOutputInnerType], Output] = at[Array[MyriadOutputInnerType]] { types =>
    (evalFunction,schemaDefRequirement) =>
      types.toList.map(_.fold(MyriadOutputInnerTypeToWomValue).apply(evalFunction, schemaDefRequirement)).map(_.toErrorOr).collectFirst({
        case Valid(validValue) => validValue
      }).toValidNel(s"Cannot find a suitable type to build a WomValue from in ${types.mkString(", ")}").toIOChecked
  }
}

object MyriadOutputInnerTypeToWomValue extends Poly1 {

  import Case._
  import MyriadOutputTypeToWomValue.Output

  def ex(component: String) = throw new RuntimeException(s"output type $component cannot yield a wom value")

  implicit def cwlType: Aux[CwlType, Output] = at[CwlType] { _ => (_,_) =>
    "No output binding is defined. Are you expecting the output to be inferred from a cwl.output.json file ? If so please make sure the file was effectively created.".invalidIOChecked
  }

  implicit def ors: Aux[OutputRecordSchema, Output] = at[OutputRecordSchema] { ors => (evalFunction, schemaDefRequirement) => ors match {
    case OutputRecordSchema(_, Some(fields), _) =>
        // Go over each field and evaluate the binding if there's one, otherwise keep folding over field types
        def evaluateValues = fields.toList.traverse[IOChecked, ((String, WomValue), (String, WomType))]({ field =>
          val womType = field.`type`.fold(MyriadOutputTypeToWomType).apply(schemaDefRequirement)
          val womValue: IOChecked[WomValue] = field.outputBinding match {
            case Some(binding) => evalFunction(binding, womType)
            case None => field.`type`.fold(MyriadOutputTypeToWomValue).apply(evalFunction, schemaDefRequirement)
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
          WomObject.withTypeUnsafe(valueMap.toMap, WomCompositeType(typeMap.toMap))
        }
    case ors => ors.toString |> ex
  }}

  implicit def oes: Aux[OutputEnumSchema, Output] = at[OutputEnumSchema]{
    //DB: I tried to do a pattern match as the overall function here but the compiler exploded
    oes => (f, schemaDefRequirement) => oes match {
      case oes@OutputEnumSchema(_, _, _, _, Some(outputBinding)) => f(outputBinding, oes.toWomEnumerationType)
      case _ => s"The enumeration type $oes requires an outputbinding to be evaluated.".invalidIOChecked
    }
  }

  implicit def oas: Aux[OutputArraySchema, Output] = at[OutputArraySchema]{ oas => (evalFunction, schemaDefRequirement) => oas match {
    case OutputArraySchema(itemsType, _, _, outputBinding) =>
      lazy val itemsWomType = itemsType.fold(MyriadOutputTypeToWomType).apply(schemaDefRequirement)
      def fromBinding = outputBinding.map(evalFunction(_, WomArrayType(itemsWomType)))
      def fromTypes = itemsType.fold(MyriadOutputTypeToWomValue).apply(evalFunction, schemaDefRequirement)
      fromBinding.getOrElse(fromTypes)
  }}

  implicit def s: Aux[String, Output] = at[String]{ s => (_, _) =>
    s.toString |> ex
  }
}
