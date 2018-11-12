package cwl

import cats.instances.list._
import cats.syntax.traverse._
import common.validation.IOChecked._
import cwl.CwlType.CwlType
import cwl.MyriadOutputTypeToWomFiles.EvaluationFunction
import mouse.all._
import shapeless.Poly1
import wom.expression.FileEvaluation

object MyriadOutputTypeToWomFiles extends Poly1 {

  type EvaluationFunction = CommandOutputBinding => IOChecked[Set[FileEvaluation]]

  import Case._

  implicit def cwlType: Aux[MyriadOutputInnerType, EvaluationFunction => IOChecked[Set[FileEvaluation]]] = at[MyriadOutputInnerType]{
    _.fold(MyriadOutputInnerTypeToWomFiles)
  }

  implicit def acwl: Aux[Array[MyriadOutputInnerType], EvaluationFunction => IOChecked[Set[FileEvaluation]]] = at[Array[MyriadOutputInnerType]] { types =>
    evalFunction =>
      types.toList.traverse(_.fold(MyriadOutputInnerTypeToWomFiles).apply(evalFunction)).map(_.toSet.flatten)
  }
}

object MyriadOutputInnerTypeToWomFiles extends Poly1 {

  import Case._

  def ex(component: String) = throw new RuntimeException(s"output type $component cannot yield wom files")

  implicit def cwlType: Aux[CwlType, EvaluationFunction => IOChecked[Set[FileEvaluation]]] = at[CwlType] { _ =>_ =>
    Set.empty[FileEvaluation].validIOChecked
  }

  implicit def ors: Aux[OutputRecordSchema, EvaluationFunction => IOChecked[Set[FileEvaluation]]] = at[OutputRecordSchema] {
    case OutputRecordSchema(_, Some(fields), _) =>
      evalFunction =>
        fields.toList.traverse[IOChecked, Set[FileEvaluation]]({ field =>
          field.outputBinding match {
            case Some(binding) => evalFunction(binding)
            case None => field.`type`.fold(MyriadOutputTypeToWomFiles).apply(evalFunction)
          }
        }).map(_.toSet.flatten)
    case ors => ors.toString |> ex
  }

  implicit def oes: Aux[OutputEnumSchema, EvaluationFunction => IOChecked[Set[FileEvaluation]]] = at[OutputEnumSchema]{ oes =>_ =>
    oes.toString |> ex
  }

  implicit def oas: Aux[OutputArraySchema, EvaluationFunction => IOChecked[Set[FileEvaluation]]] = at[OutputArraySchema]{ oas =>
    evalFunction =>
      import cats.syntax.apply._
      def fromBinding: IOChecked[Set[FileEvaluation]] = oas.outputBinding.map(evalFunction).getOrElse(Set.empty[FileEvaluation].validIOChecked)
      def fromType: IOChecked[Set[FileEvaluation]] = oas.items.fold(MyriadOutputTypeToWomFiles).apply(evalFunction)

      (fromBinding, fromType) mapN (_ ++ _)
  }

  implicit def s: Aux[String, EvaluationFunction => IOChecked[Set[FileEvaluation]]] = at[String]{ _ =>_ =>
    Set.empty[FileEvaluation].validIOChecked
  }
}
