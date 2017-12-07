package cwl

import cwl.CwlType.CwlType
import shapeless.{Coproduct, Poly1}
import wom.types.{WomArrayType, WomType}
import cats.syntax.foldable._
import cats.instances.list._
import cats.instances.set._

object MyriadOutputTypeToWomType extends Poly1{

  import Case._

  implicit def cwlType: Aux[CwlType, WomType] = at[CwlType]{
    Coproduct[MyriadOutputInnerType](_).fold(MyriadOutputArrayTypeToWomType)
  }

  implicit def ors: Aux[OutputRecordSchema, WomType] = at[OutputRecordSchema] {
    Coproduct[MyriadOutputInnerType](_).fold(MyriadOutputArrayTypeToWomType)
  }

  implicit def oes: Aux[OutputEnumSchema, WomType] = at[OutputEnumSchema] {
    Coproduct[MyriadOutputInnerType](_).fold(MyriadOutputArrayTypeToWomType)
  }

  implicit def oas: Aux[OutputArraySchema, WomType] = at[OutputArraySchema] {
    Coproduct[MyriadOutputInnerType](_).fold(MyriadOutputArrayTypeToWomType)
  }

  implicit def s: Aux[String, WomType] =  at[String] {
    Coproduct[MyriadOutputInnerType](_).fold(MyriadOutputArrayTypeToWomType)
  }

  implicit def acwl: Aux[Array[MyriadOutputInnerType], WomType] = at[Array[MyriadOutputInnerType]] {
    _.toList.foldMap(i => Set(i)).map(_.fold(MyriadOutputArrayTypeToWomType)).toList match {
      case head :: Nil => WomArrayType(head)
      case _ => throw new RuntimeException("Wom does not provide an array of >1 types")
    }
  }
}

object MyriadOutputArrayTypeToWomType extends Poly1 {

  import Case._

  implicit def cwlType: Aux[CwlType, WomType] = at[CwlType]{
    cwl.cwlTypeToWdlType
  }

  implicit def ors: Aux[OutputRecordSchema, WomType] = at[OutputRecordSchema] {
    _ => ???
  }

  implicit def oes: Aux[OutputEnumSchema, WomType] = at[OutputEnumSchema] {_ => ???}

  implicit def oas: Aux[OutputArraySchema, WomType] = at[OutputArraySchema] {
    oas =>
      val arrayType: WomType = oas.items.fold(MyriadOutputTypeToWomType)

      WomArrayType(arrayType)
  }
  implicit def s: Aux[String, WomType] =  at[String] { _ => ???}
}
