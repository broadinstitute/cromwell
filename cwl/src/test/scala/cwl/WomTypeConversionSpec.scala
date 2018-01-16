package cwl

import org.scalacheck.Properties
import org.scalacheck.Prop._
import shapeless.Coproduct
import wom.types.{WomArrayType, WomOptionalType, WomStringType}
import mouse.`try`._

import scala.util.Try

class WomTypeConversionSpec extends Properties("CWL -> WOM Conversion"){

  property("ArrayInputSchema") = secure {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    val x = InputArraySchema(items = Coproduct[MyriadInputType](y))
    val z = Coproduct[MyriadInputInnerType](x)
    Coproduct[MyriadInputType](z).fold(MyriadInputTypeToWomType) == WomArrayType(WomStringType)
  }

  property("Cwl String") = secure {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    Coproduct[MyriadInputType](y).fold(MyriadInputTypeToWomType) == WomStringType
  }

  property("Array of a single type is actually one type not in an array") = secure {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    Coproduct[MyriadInputType](Array(y)).fold(MyriadInputTypeToWomType) == WomStringType
  }

  property("Array of more than one type is not allowed") = secure {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    val z = Coproduct[MyriadInputInnerType](CwlType.Boolean)
    Try(Coproduct[MyriadInputType](Array(y, z)).fold(MyriadInputTypeToWomType)).cata(
      s =>  throw new RuntimeException(s"failure expected but got $s!"),
      _.getMessage == "Multi types not supported yet"
    )
  }

  property("Array of a type followed by a null is an optional type") = secure {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    val z = Coproduct[MyriadInputInnerType](CwlType.Null)
    Try(Coproduct[MyriadInputType](Array(y, z)).fold(MyriadInputTypeToWomType)).cata(
      success => (success == WomOptionalType(WomStringType)) :| "input should evaluate to an optional string type",
      failure => false :| s"expected success but received a failure"
    )
  }

  property("ArrayOutputSchema") = secure {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    val x = OutputArraySchema(items = Coproduct[MyriadOutputType](y))
    val z = Coproduct[MyriadOutputInnerType](x)
    Coproduct[MyriadOutputType](z).fold(MyriadOutputTypeToWomType) == WomArrayType(WomStringType)
  }

  property("Cwl String") = secure {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    Coproduct[MyriadOutputType](y).fold(MyriadOutputTypeToWomType) == WomStringType
  }

  property("Array of a single type is the same as one type not in an array") = secure {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    Coproduct[MyriadOutputType](Array(y)).fold(MyriadOutputTypeToWomType) == WomStringType
  }

  property("Array of more than one type is not allowed") = throws(classOf[NotImplementedError]) {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    val z = Coproduct[MyriadOutputInnerType](CwlType.Boolean)
    Coproduct[MyriadOutputType](Array(y, z)).fold(MyriadOutputTypeToWomType)
  }

  property("Array of a type followed by a null is an optional type") = secure {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    val z = Coproduct[MyriadOutputInnerType](CwlType.Null)
    Try(Coproduct[MyriadOutputType](Array(y, z)).fold(MyriadOutputTypeToWomType)).cata(
      success => (success == WomOptionalType(WomStringType)) :| "input should evaluate to an optional string type",
      failure => false :| s"expected success but received a failure"
    )
  }
}
