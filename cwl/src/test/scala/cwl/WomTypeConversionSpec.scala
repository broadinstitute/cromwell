package cwl

import org.scalacheck.Properties
import org.scalacheck.Prop._
import shapeless.Coproduct
import wom.types.{WomArrayType, WomStringType}

class WomTypeConversionSpec extends Properties("CWL -> WOM Conversion"){

  property("ArrayInputSchema") = secure {
    val y = Coproduct[MyriadInputInnerType ](CwlType.String)
    val x = InputArraySchema(items = Coproduct[MyriadInputType](y))
    val z = Coproduct[MyriadInputInnerType](x)
    Coproduct[MyriadInputType](z).fold(MyriadInputTypeToWomType) == WomArrayType(WomStringType)
  }

  property("Cwl String") = secure {
    val y = Coproduct[MyriadInputInnerType ](CwlType.String)
    Coproduct[MyriadInputType](y).fold(MyriadInputTypeToWomType) == WomStringType
  }

  property("Array of things is array of one type") = secure {
    val y = Coproduct[MyriadInputInnerType ](CwlType.String)
    Coproduct[MyriadInputType](Array(y)).fold(MyriadInputTypeToWomType) == WomArrayType(WomStringType)
  }

  property("Array of more than one type is not allowed") = throws(classOf[RuntimeException]) {
    val y = Coproduct[MyriadInputInnerType ](CwlType.String)
    val z = Coproduct[MyriadInputInnerType ](CwlType.Boolean)
    Coproduct[MyriadInputType](Array(y, z)).fold(MyriadInputTypeToWomType)
  }
}
