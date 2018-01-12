package cwl

import cwl.CwlType.CwlType
import cwl.command.ParentName
import mouse.all._
import shapeless.Poly1
import wom.types._


object MyriadInputTypeToWomType extends Poly1 {

  implicit def m = at[MyriadInputInnerType] {_.fold(MyriadInputInnerTypeToWomType)}

  // An array of type means "this input value can be in any of those types"
  implicit def am = at[Array[MyriadInputInnerType]] { types =>
    types.partition(_.select[CwlType].contains(CwlType.Null)) match {
        // If there's a single non null type, use that
      case (nullTypes, Array(singleNonNullType)) if nullTypes.isEmpty =>
        singleNonNullType.fold(MyriadInputInnerTypeToWomType)
        // If there's a null type and a single non null type, it's a WomOptionalType
      case (nullTypes, Array(singleNonNullType)) if nullTypes.nonEmpty =>
        WomOptionalType(singleNonNullType.fold(MyriadInputInnerTypeToWomType))
        // Leave other "Coproduct types" unsupported for now
      case _ =>
        throw new NotImplementedError("Multi types not supported yet")
    }
  }
}

object MyriadInputInnerTypeToWomType extends Poly1 {
  import Case._

  def ex(component: String) = throw new RuntimeException(s"input type $component not yet suported by WOM!")

  implicit def ct: Aux[CwlType, WomType] = at[CwlType]{
    cwl.cwlTypeToWomType
  }
  implicit def irs: Aux[InputRecordSchema, WomType] = at[InputRecordSchema]{
    case InputRecordSchema(_, Some(fields), _) =>
      val typeMap = fields.map({ field =>
          FullyQualifiedName(field.name)(ParentName.empty).id -> field.`type`.fold(MyriadInputTypeToWomType)
      }).toMap
      WomCompositeType(typeMap)
    case irs => irs.toString |> ex
  }

  implicit def ies: Aux[InputEnumSchema, WomType] = at[InputEnumSchema]{
    _.toString |> ex
  }

  implicit def ias: Aux[InputArraySchema, WomType] = at[InputArraySchema]{
    ias =>
      val arrayType: WomType = ias.items.fold(MyriadInputTypeToWomType)

      WomArrayType(arrayType)
  }
  implicit def s: Aux[String, WomType] = at[String]{
    _.toString |> ex
  }

}
