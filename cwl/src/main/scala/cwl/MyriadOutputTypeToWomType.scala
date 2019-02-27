package cwl

import cats.data.NonEmptyList
import cwl.CwlType.CwlType
import cwl.command.ParentName
import mouse.all._
import shapeless.Poly1
import wom.types._

object MyriadOutputInnerTypeToString extends Poly1 {

  implicit def cwlType = at[CwlType] {
    _.toString
  }

  implicit def ors = at[OutputRecordSchema] {
    _.toString
  }

  implicit def oes = at[OutputEnumSchema] {
    _.toString
  }

  implicit def oas = at[OutputArraySchema] {
    _.toString
  }

  implicit def s = at[String] {
    identity
  }
}

object MyriadOutputTypeToWomType extends Poly1{

  import Case._

  type SchemaDefToWomType = SchemaDefRequirement => WomType

  implicit def cwlType: Aux[MyriadOutputInnerType, SchemaDefToWomType] = at[MyriadOutputInnerType]{ moit => schemaDefRequirement =>
    moit.fold(MyriadOutputInnerTypeToWomType).apply(schemaDefRequirement)
  }

  implicit def acwl: Aux[Array[MyriadOutputInnerType], SchemaDefToWomType] = at[Array[MyriadOutputInnerType]] { types => schemaDefRequirement =>
    types.partition(_.select[CwlType].contains(CwlType.Null)) match {
      // If there's a single non null type, use that
      case (Array(), Array(singleNonNullType)) =>
        singleNonNullType.fold(MyriadOutputInnerTypeToWomType).apply(schemaDefRequirement)
      case (Array(), array: Array[MyriadOutputInnerType]) if array.size > 1  =>
        val types = array.map(_.fold(MyriadOutputInnerTypeToWomType).apply(schemaDefRequirement))
        WomCoproductType(NonEmptyList.fromListUnsafe(types.toList))
      // If there's a null type and a single non null type, it's a WomOptionalType
      case (Array(_), Array(singleNonNullType)) =>
        WomOptionalType(singleNonNullType.fold(MyriadOutputInnerTypeToWomType).apply(schemaDefRequirement))
      case (Array(_), array: Array[MyriadOutputInnerType]) if array.size > 1  =>
        val types = array.map(_.fold(MyriadOutputInnerTypeToWomType).apply(schemaDefRequirement))
        WomOptionalType(WomCoproductType(NonEmptyList.fromListUnsafe(types.toList)))
      case _ =>
        val readableTypes = types.map(_.fold(MyriadOutputInnerTypeToString)).mkString(", ")
        throw new UnsupportedOperationException(s"Cromwell only supports single types or optionals (as indicated by [null, X]). Instead we saw: $readableTypes")
    }
  }
}

object MyriadOutputInnerTypeToWomType extends Poly1 {

  import Case._
  import MyriadOutputTypeToWomType.SchemaDefToWomType

  implicit def cwlType: Aux[CwlType, SchemaDefToWomType] =
    at[CwlType]{
      cwl.cwlTypeToWomType andThen Function.const
    }

  implicit def ors: Aux[OutputRecordSchema, SchemaDefToWomType] = at[OutputRecordSchema] {
    ors => schemaDefRequirement =>
      ors.fields.fold(
        WomCompositeType(Map.empty))(
        {fields =>
          val typeMap = fields.map({ field =>
            val parsedName = FullyQualifiedName(field.name)(ParentName.empty).id
            parsedName -> field.`type`.fold(MyriadOutputTypeToWomType).apply(schemaDefRequirement)
          }).toMap
          WomCompositeType(typeMap)
        }
      )
  }

  implicit def oes: Aux[OutputEnumSchema, MyriadOutputTypeToWomType.SchemaDefToWomType] =
    at[OutputEnumSchema] {
      oes => oes.toWomEnumerationType |> Function.const
    }

  implicit def oas: Aux[OutputArraySchema, MyriadOutputTypeToWomType.SchemaDefToWomType] =
    at[OutputArraySchema] {
      oas => schemaDefRequirement =>
        val arrayType: WomType = oas.items.fold(MyriadOutputTypeToWomType).apply(schemaDefRequirement)
        WomArrayType(arrayType)
    }

  implicit def s: Aux[String, MyriadOutputTypeToWomType.SchemaDefToWomType] =
    at[String] {
      string => schemaReq =>
        schemaReq.lookupType(string).getOrElse(throw new RuntimeException(s"Custom type $string was referred to but not found in schema def ${schemaReq.types.mkString(", ")}."))
    }
}
