package cwl

import cwl.CwlType.CwlType
import mouse.all._
import shapeless.Poly1

object MyriadOutputTypeToGlob extends Poly1 {

  import Case._

  implicit def cwlType: Aux[MyriadOutputInnerType, List[Glob]] = at[MyriadOutputInnerType]{
    _.fold(MyriadOutputInnerTypeToGlob)
  }

  implicit def acwl: Aux[Array[MyriadOutputInnerType], List[Glob]] = at[Array[MyriadOutputInnerType]] { 
    _.toList.flatMap(_.fold(MyriadOutputInnerTypeToGlob))
  }
}

object MyriadOutputInnerTypeToGlob extends Poly1 {

  import Case._

  def ex(component: String) = throw new RuntimeException(s"output type $component cannot yield a wom value")

  implicit def cwlType: Aux[CwlType, List[Glob]] = at[CwlType] { 
    Function.const(List.empty)
  }

  implicit def ors: Aux[OutputRecordSchema, List[Glob]] = at[OutputRecordSchema] {
    case OutputRecordSchema(_, Some(fields), _) =>
        fields.toList.flatMap({ field =>
          field.outputBinding.flatMap(_.glob).toList ++ field.`type`.fold(MyriadOutputTypeToGlob)
        })
    case ors => ors.toString |> ex
  }

  implicit def oes: Aux[OutputEnumSchema, List[Glob]] = at[OutputEnumSchema]{ oes =>
    oes.toString |> ex
  }

  implicit def oas: Aux[OutputArraySchema, List[Glob]] = at[OutputArraySchema]{ oas =>
    oas.outputBinding.flatMap(_.glob).toList ++ oas.items.fold(MyriadOutputTypeToGlob)
  }

  implicit def s: Aux[String, List[Glob]] = at[String]{
    Function.const(List.empty)
  }
}
