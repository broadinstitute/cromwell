package cwl

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import cwl.CwlType.CwlType
import cwl.MyriadOutputTypeToWomFiles.EvaluationFunction
import mouse.all._
import shapeless.Poly1
import wom.values.WomFile

object MyriadOutputTypeToWomFiles extends Poly1 {

  type EvaluationFunction = CommandOutputBinding => ErrorOr[Set[WomFile]]

  import Case._

  implicit def cwlType: Aux[MyriadOutputInnerType, EvaluationFunction => ErrorOr[Set[WomFile]]] = at[MyriadOutputInnerType]{
    _.fold(MyriadOutputInnerTypeToWomFiles)
  }

  implicit def acwl: Aux[Array[MyriadOutputInnerType], EvaluationFunction => ErrorOr[Set[WomFile]]] = at[Array[MyriadOutputInnerType]] { types => 
    evalFunction =>
      types.toList.traverse[ErrorOr, Set[WomFile]](_.fold(MyriadOutputInnerTypeToWomFiles).apply(evalFunction)).map(_.toSet.flatten)
  }
}

object MyriadOutputInnerTypeToWomFiles extends Poly1 {

  import Case._

  def ex(component: String) = throw new RuntimeException(s"output type $component cannot yield wom files")

  implicit def cwlType: Aux[CwlType, EvaluationFunction => ErrorOr[Set[WomFile]]] = at[CwlType] { _ => _ =>
    Set.empty[WomFile].validNel
  }

  implicit def ors: Aux[OutputRecordSchema, EvaluationFunction => ErrorOr[Set[WomFile]]] = at[OutputRecordSchema] {
    case OutputRecordSchema(_, Some(fields), _) =>
      evalFunction =>
        fields.toList.traverse[ErrorOr, Set[WomFile]]({ field =>
          field.outputBinding match {
            case Some(binding) => evalFunction(binding)
            case None => field.`type`.fold(MyriadOutputTypeToWomFiles).apply(evalFunction)
          }
        }).map(_.toSet.flatten)
    case ors => ors.toString |> ex
  }

  implicit def oes: Aux[OutputEnumSchema, EvaluationFunction => ErrorOr[Set[WomFile]]] = at[OutputEnumSchema]{ oes => _ =>
    oes.toString |> ex
  }

  implicit def oas: Aux[OutputArraySchema, EvaluationFunction => ErrorOr[Set[WomFile]]] = at[OutputArraySchema]{ oas =>
    evalFunction =>
      import cats.syntax.apply._
      def fromBinding: ErrorOr[Set[WomFile]] = oas.outputBinding.map(evalFunction).getOrElse(Set.empty[WomFile].validNel)
      def fromType: ErrorOr[Set[WomFile]] = oas.items.fold(MyriadOutputTypeToWomFiles).apply(evalFunction)

      (fromBinding, fromType) mapN (_ ++ _)
  }

  implicit def s: Aux[String, EvaluationFunction => ErrorOr[Set[WomFile]]] = at[String]{ _ => _ =>
    Set.empty[WomFile].validNel
  }
}
