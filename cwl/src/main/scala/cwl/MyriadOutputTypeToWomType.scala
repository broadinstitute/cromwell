package cwl

import cwl.CwlType.CwlType
import shapeless.Poly1
import wom.types.{WomArrayType, WomType}
import cats.syntax.foldable._
import cats.instances.list._
import cats.instances.set._
import mouse.all._

object MyriadOutputTypeToWomType extends Poly1{

  import Case._

  implicit def cwlType: Aux[MyriadOutputInnerType, WomType] = at[MyriadOutputInnerType]{
    _.fold(MyriadOutputInnerTypeToWomType)
  }

  implicit def acwl: Aux[Array[MyriadOutputInnerType], WomType] = at[Array[MyriadOutputInnerType]] {
    _.toList.foldMap(i => Set(i)).map(_.fold(MyriadOutputInnerTypeToWomType)).toList match {
      case head :: Nil => WomArrayType(head)
      case _ => throw new RuntimeException("Wom does not provide an array of >1 types")
    }
  }
}

object MyriadOutputInnerTypeToWomType extends Poly1 {

  import Case._

  def ex(component: String) = throw new RuntimeException(s"input type $component not yet suported by WOM!")

  implicit def cwlType: Aux[CwlType, WomType] = at[CwlType]{
    cwl.cwlTypeToWomType
  }

  implicit def ors: Aux[OutputRecordSchema, WomType] = at[OutputRecordSchema] {
    _.toString |> ex
  }

  implicit def oes: Aux[OutputEnumSchema, WomType] = at[OutputEnumSchema] {
    _.toString |> ex
  }

  implicit def oas: Aux[OutputArraySchema, WomType] = at[OutputArraySchema] {
    oas =>
      val arrayType: WomType = oas.items.fold(MyriadOutputTypeToWomType)

      WomArrayType(arrayType)
  }
  implicit def s: Aux[String, WomType] =  at[String] { 
    _.toString |> ex
  }
}
