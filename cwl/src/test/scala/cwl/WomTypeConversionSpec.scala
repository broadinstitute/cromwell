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

  property("Array of more than one type is not allowed") = throws(classOf[NotImplementedError]){
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    val z = Coproduct[MyriadInputInnerType](CwlType.Boolean)
    Coproduct[MyriadInputType](Array(y, z)).fold(MyriadInputTypeToWomType)
  }

  property("a 2-element Array of a single type accompanied by a null is an optional type") = secure {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    val z = Coproduct[MyriadInputInnerType](CwlType.Null)
    testInputArray(Array(y,z), WomOptionalType(WomStringType))
  }


  def testInputArray(array: Array[MyriadInputInnerType], assertedType: WomType ) = {
    def f(array: Array[MyriadInputInnerType]) =
      Try(Coproduct[MyriadInputType](array).fold(MyriadInputTypeToWomType)).cata(
        success => (success == assertedType) :| "input should evaluate to an optional array of files type",
        failure => false :| s"expected an optional array of files type but received a failure"
      )
    f(array) && f(array.reverse)
  }

  property("Optional Array of a type is interpreted correctly as input type") = secure {
    val miit = Coproduct[MyriadInputInnerType](CwlType.File)
    val mit = Coproduct[MyriadInputType](miit)
    val ias =  InputArraySchema(items = mit)
    val y = Coproduct[MyriadInputInnerType](ias)
    val z = Coproduct[MyriadInputInnerType](CwlType.Null)
    testInputArray(Array(y, z), WomOptionalType(WomMaybeEmptyArrayType(WomMaybePopulatedFileType)))
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

  property("Output Array of more than one type is not allowed") = throws(classOf[NotImplementedError]) {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    val z = Coproduct[MyriadOutputInnerType](CwlType.Boolean)
    Coproduct[MyriadOutputType](Array(y, z)).fold(MyriadOutputTypeToWomType)
  }

  def testOutputArray(array: Array[MyriadOutputInnerType], assertedType: WomType) = {
    def f(array: Array[MyriadOutputInnerType]) =
      Try(Coproduct[MyriadOutputType](array).fold(MyriadOutputTypeToWomType)).cata(
        success => (success == assertedType) :| "input should evaluate to an optional string type",
        failure => false :| s"expected an optional string type but received a failure"
      )
    f(array) && f(array.reverse)
  }

  property("a 2-element Array of a single type accompanied by a null is an optional type") = secure {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    val z = Coproduct[MyriadOutputInnerType](CwlType.Null)
    testOutputArray(Array(y,z), WomOptionalType(WomStringType))
  }

  property("Optional Array of a type is interpreted correctly as output type") = secure {
    val moit = Coproduct[MyriadOutputInnerType](CwlType.File)
    val mot = Coproduct[MyriadOutputType](moit)
    val oas =  OutputArraySchema(items = mot)
    val moit2 = Coproduct[MyriadOutputInnerType](oas)
    val moit3 = Coproduct[MyriadOutputInnerType](CwlType.Null)
    testOutputArray(Array(moit2,moit3), WomOptionalType(WomMaybeEmptyArrayType(WomMaybePopulatedFileType)))
  }
}
