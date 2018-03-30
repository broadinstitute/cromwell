package cwl

import cats.data.NonEmptyList
import cwl.CwlType.CwlType
import cwl.MyriadInputTypeToWomType.SchemaLookup
import cwl.command.ParentName
import mouse.all._
import shapeless.{Inl, Inr, Poly1}
import wom.types._

object MyriadInputInnerTypeToString extends Poly1 {
  implicit def ct = at[CwlType]{ _.toString }
  implicit def irs = at[InputRecordSchema]{_.toString}
  implicit def ies = at[InputEnumSchema]{ _.toString }
  implicit def ias = at[InputArraySchema]{ _.toString}
  implicit def s = at[String]{identity}
}

object MyriadInputTypeToWomType extends Poly1 {

  import Case._

  type SchemaLookup = SchemaDefRequirement => WomType

  implicit def m:Aux[MyriadInputInnerType, SchemaLookup]= at[MyriadInputInnerType] {_.fold(MyriadInputInnerTypeToWomType)}

  // An array of type means "this input value can be in any of those types."
  // Currently we only accept single types or [null, X] to mean Optional[X]
  implicit def am: Aux[Array[MyriadInputInnerType], SchemaLookup] = at[Array[MyriadInputInnerType]] {
    types =>
      schemaLookup =>
        types.partition(_.select[CwlType].contains(CwlType.Null)) match {
          // If there's a single non null type, use that
          case (Array(), Array(singleNonNullType)) =>
            singleNonNullType.fold(MyriadInputInnerTypeToWomType).apply(schemaLookup)
          case (Array(), array: Array[MyriadInputInnerType]) if array.length > 1  =>
            val types = array.map(_.fold(MyriadInputInnerTypeToWomType).apply(schemaLookup))
            WomCoproductType(NonEmptyList.fromListUnsafe(types.toList))
          // If there's a null type and a single non null type, it's a WomOptionalType
          case (Array(_), Array(singleNonNullType)) =>
            WomOptionalType(singleNonNullType.fold(MyriadInputInnerTypeToWomType).apply(schemaLookup))
          case (Array(_), array: Array[MyriadInputInnerType]) if array.length > 1  =>
            val types = array.map(_.fold(MyriadInputInnerTypeToWomType).apply(schemaLookup))
            WomOptionalType(WomCoproductType(NonEmptyList.fromListUnsafe(types.toList)))
          case (Array(_), _) =>
            val readableTypes = types.map(_.fold(MyriadInputInnerTypeToString)).mkString(", ")
            throw new NotImplementedError(s"Cromwell only supports single types or optionals (as indicated by [null, X]). Instead we saw: $readableTypes")
    }
  }
}

object MyriadInputInnerTypeToWomType extends Poly1 {
  import Case._

  def ex(component: String) = throw new RuntimeException(s"input type $component not yet suported by WOM!")

  implicit def ct: Aux[CwlType, SchemaLookup] = at[CwlType]{
    ct => _ =>
    cwl.cwlTypeToWomType(ct)
  }

  def inputRecordSchemaToWomType(irs: InputRecordSchema): SchemaLookup = { schemaLookup: SchemaDefRequirement =>
    irs match {
      case InputRecordSchema(_, _, Some(fields), _) =>
        val typeMap = fields.map({ field =>
          FullyQualifiedName(field.name)(ParentName.empty).id -> field.`type`.fold(MyriadInputTypeToWomType).apply(schemaLookup)
        }).toMap
        WomCompositeType(typeMap)
      case irs => irs.toString |> ex
    }
  }

  implicit def irs: Aux[InputRecordSchema, SchemaLookup] = at[InputRecordSchema]{
    inputRecordSchemaToWomType
  }

  implicit def ies: Aux[InputEnumSchema, SchemaLookup] = at[InputEnumSchema]{
    ies => _ =>
      ies.toWomEnumerationType
  }

  implicit def ias: Aux[InputArraySchema, SchemaLookup] = at[InputArraySchema]{
    ias =>
      lookup =>
        val arrayType: WomType = ias.items.fold(MyriadInputTypeToWomType).apply(lookup)
        WomArrayType(arrayType)
  }
  implicit def s: Aux[String, SchemaLookup] = at[String]{
    string =>
      schemaReq =>
        schemaReq.types.toList.flatMap{
          case Inl(inputRecordSchema: InputRecordSchema) if inputRecordSchema.name == string.drop(1) =>
            List(inputRecordSchemaToWomType(inputRecordSchema).apply(schemaReq))
          case _ => List()
        }.headOption.getOrElse(throw new RuntimeException(s"Custom type $string was referred to but not found."))
  }

}
