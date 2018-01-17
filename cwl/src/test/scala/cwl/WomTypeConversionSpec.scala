package cwl

import org.scalacheck.Properties
import org.scalacheck.Prop._
import shapeless.Coproduct
import wom.types._
import mouse.`try`._

import scala.util.Try

class WomTypeConversionSpec extends Properties("CWL -> WOM Conversion"){

  /* ******** Inputs *********** */
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
      failure => false :| s"expected optional string type but received a failure"
    )
  }

  property("Optional Array of a type is interpreted correctly as input type") = secure {
    val miit = Coproduct[MyriadInputInnerType](CwlType.File)
    val mit = Coproduct[MyriadInputType](miit)
    val ias =  InputArraySchema(items = mit)
    val y = Coproduct[MyriadInputInnerType](ias)
    val z = Coproduct[MyriadInputInnerType](CwlType.Null)
    def testInputArray(array: Array[MyriadInputInnerType]) = {
      Try(Coproduct[MyriadInputType](array).fold(MyriadInputTypeToWomType)).cata(
        success => (success == WomOptionalType(WomMaybeEmptyArrayType(WomMaybePopulatedFileType))) :| "input should evaluate to an optional array of files type",
        failure => false :| s"expected an optional array of files type but received a failure"
      )
    }
    testInputArray(Array(y, z)) && testInputArray(Array(z,y))
  }

  /* ******** Outputs *********** */
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
      failure => false :| s"expected an optional string type but received a failure"
    )
  }

  property("Optional Array of a type is interpreted correctly as output type") = secure {
    def testOutputArray(array: Array[MyriadOutputInnerType]) = {
      Try(Coproduct[MyriadOutputType](array).fold(MyriadOutputTypeToWomType)).cata(
        success => (success == WomOptionalType(WomMaybeEmptyArrayType(WomMaybePopulatedFileType))) :| "input should evaluate to an optional array of files type",
        failure => false :| s"expected to receive an optional array of files type but received a failure"
      )
    }
    val miit = Coproduct[MyriadOutputInnerType](CwlType.File)
    val mit = Coproduct[MyriadOutputType](miit)
    val ias =  OutputArraySchema(items = mit)
    val y = Coproduct[MyriadOutputInnerType](ias)
    val z = Coproduct[MyriadOutputInnerType](CwlType.Null)
    testOutputArray(Array(y,z)) && testOutputArray(Array(z,y))
  }
}
