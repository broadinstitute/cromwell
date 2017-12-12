package cwl

import org.scalacheck.Properties
import org.scalacheck.Prop._
import shapeless.Coproduct
import wom.types.{WomArrayType, WomStringType}

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

  property("Array of things is array of one type") = secure {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    Coproduct[MyriadInputType](Array(y)).fold(MyriadInputTypeToWomType) == WomArrayType(WomStringType)
  }

  property("Array of more than one type is not allowed") = secure {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    val z = Coproduct[MyriadInputInnerType](CwlType.Boolean)
    Try(Coproduct[MyriadInputType](Array(y, z)).fold(MyriadInputTypeToWomType)).cata(
      s =>  throw new RuntimeException(s"failure expected but got $s!"),
      _.getMessage == "Wom does not provide an array of >1 types"
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

  property("Array of things is array of one type") = secure {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    Coproduct[MyriadOutputType](Array(y)).fold(MyriadOutputTypeToWomType) == WomArrayType(WomStringType)
  }

  property("Array of more than one type is not allowed") = throws(classOf[RuntimeException]) {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    val z = Coproduct[MyriadOutputInnerType](CwlType.Boolean)
    Coproduct[MyriadOutputType](Array(y, z)).fold(MyriadOutputTypeToWomType)
  }
}
