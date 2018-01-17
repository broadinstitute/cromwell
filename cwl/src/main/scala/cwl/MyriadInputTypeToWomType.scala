package cwl

import cwl.CwlType.CwlType
import cwl.command.ParentName
import mouse.all._
import shapeless.Poly1
import wom.types._

object MyriadInputInnerTypeToString extends Poly1 {
  implicit def ct = at[CwlType]{ _.toString }
  implicit def irs = at[InputRecordSchema]{_.toString}
  implicit def ies = at[InputEnumSchema]{ _.toString }
  implicit def ias = at[InputArraySchema]{ _.toString}
  implicit def s = at[String]{identity}
}

object MyriadInputTypeToWomType extends Poly1 {

  implicit def m = at[MyriadInputInnerType] {_.fold(MyriadInputInnerTypeToWomType)}

  // An array of type means "this input value can be in any of those types."
  // Currently we only accept single types or [null, X] to mean Optional[X]
  implicit def am = at[Array[MyriadInputInnerType]] { types =>
    types.partition(_.select[CwlType].contains(CwlType.Null)) match {
        // If there's a single non null type, use that
      case (Array(), Array(singleNonNullType)) =>
        singleNonNullType.fold(MyriadInputInnerTypeToWomType)
        // If there's a null type and a single non null type, it's a WomOptionalType
      case (Array(_), Array(singleNonNullType)) =>
        WomOptionalType(singleNonNullType.fold(MyriadInputInnerTypeToWomType))
        // Leave other "Coproduct types" unsupported for now
      case _ =>
        val readableTypes = types.map(_.fold(MyriadInputInnerTypeToString)).mkString(", ")
        throw new NotImplementedError(s"Cromwell only supports single types or optionals (as indicated by [null, X]). Instead we saw: $readableTypes")
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
