package cwl

import shapeless.Poly1
import wom.types.{WomArrayType, WomCompositeType, WomType}
import cwl.CwlType.CwlType
import cats.syntax.foldable._
import cats.instances.list._
import cats.instances.set._
import cwl.command.ParentName
import mouse.all._


object MyriadInputTypeToWomType extends Poly1 {

  implicit def m = at[MyriadInputInnerType] {_.fold(MyriadInputInnerTypeToWomType)}

  //We only accept arrays of a single type, so we have to check whether the array boils down to an array of a single type
  implicit def am = at[Array[MyriadInputInnerType]] {
    //reduce the Array down to a set to remove dupes.
    _.toList.foldMap(
      i => Set(i)
    ).map(
      _.fold(MyriadInputInnerTypeToWomType)
    ).toList match {
      case head :: Nil => WomArrayType(head)
      case _ => throw new RuntimeException("Wom does not provide an array of >1 types")
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
