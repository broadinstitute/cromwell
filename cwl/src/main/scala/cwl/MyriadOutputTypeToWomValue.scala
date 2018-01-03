package cwl

import cats.instances.list._
import cats.syntax.traverse._
import common.validation.ErrorOr.ErrorOr
import cwl.CwlType.CwlType
import cwl.MyriadOutputTypeToWomValue.EvaluationFunction
import cwl.command.ParentName
import mouse.all._
import shapeless.{Coproduct, Poly1}
import wom.types._
import wom.values.{WomArray, WomObject, WomValue}

object MyriadOutputTypeToWomValue extends Poly1 {

  type EvaluationFunction = (CommandOutputBinding, WomType) => ErrorOr[WomValue]

  import Case._

  implicit def cwlType: Aux[MyriadOutputInnerType, EvaluationFunction => ErrorOr[WomValue]] = at[MyriadOutputInnerType]{
    _.fold(MyriadOutputInnerTypeToWomValue)
  }

  implicit def acwl: Aux[Array[MyriadOutputInnerType], EvaluationFunction => ErrorOr[WomValue]] = at[Array[MyriadOutputInnerType]] { types =>
    evalFunction =>
      types.toList.traverse[ErrorOr, WomValue](_.fold(MyriadOutputInnerTypeToWomValue).apply(evalFunction)) map { values =>
        Coproduct[MyriadOutputType](types).fold(MyriadOutputTypeToWomType) match {
          case arrayType: WomArrayType => WomArray(arrayType, values)
          case other => throw new RuntimeException(s"output type $other is not an array type")
        }
      }
  }
}

object MyriadOutputInnerTypeToWomValue extends Poly1 {

  import Case._

  def ex(component: String) = throw new RuntimeException(s"output type $component cannot yield a wom value")

  implicit def cwlType: Aux[CwlType, EvaluationFunction => ErrorOr[WomValue]] = at[CwlType] { cwlType => _ =>
    cwlType.toString |> ex
  }

  implicit def ors: Aux[OutputRecordSchema, EvaluationFunction => ErrorOr[WomValue]] = at[OutputRecordSchema] {
    case OutputRecordSchema(_, Some(fields), _) =>
      evalFunction =>

        def evaluateValues = fields.toList.traverse[ErrorOr, ((String, WomValue), (String, WomType))]({ field =>
          val womType = field.`type`.fold(MyriadOutputTypeToWomType)
          val womValue: ErrorOr[WomValue] = field.outputBinding match {
            case Some(binding) =>
              evalFunction(binding, womType)
            case None => field.`type`.fold(MyriadOutputTypeToWomValue).apply(evalFunction)
          }

          womValue map { value =>
            val parsedName = FullyQualifiedName(field.name)(ParentName.empty).id
            (parsedName -> value) -> (parsedName -> womType)
          }
        })

        evaluateValues map { evaluatedValues =>
          val (valueMap, typeMap) = evaluatedValues.unzip
          WomObject.withType(valueMap.toMap, WomCompositeType(typeMap.toMap))
        }
    case ors => ors.toString |> ex
  }

  implicit def oes: Aux[OutputEnumSchema, EvaluationFunction => ErrorOr[WomValue]] = at[OutputEnumSchema]{ oes => _ =>
    oes.toString |> ex
  }

  implicit def oas: Aux[OutputArraySchema, EvaluationFunction => ErrorOr[WomValue]] = at[OutputArraySchema]{ oas => _ =>
    oas.toString |> ex
  }

  implicit def s: Aux[String, EvaluationFunction => ErrorOr[WomValue]] = at[String]{ s => _ =>
    s.toString |> ex
  }
}
