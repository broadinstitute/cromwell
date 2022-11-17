package cwl

import cats.data.NonEmptyList
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import shapeless.Coproduct
import wom.types._
import mouse.`try`._

import scala.util.Try

class WomTypeConversionSpec extends AnyFlatSpec with Matchers {
  behavior of "WomTypeConversion"

  /* ******** Inputs *********** */
  it should "convert ArrayInputSchema" in {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    val x = InputArraySchema(items = Coproduct[MyriadInputType](y))
    val z = Coproduct[MyriadInputInnerType](x)
    Coproduct[MyriadInputType](z).fold(MyriadInputTypeToWomType).apply(SchemaDefRequirement()) shouldBe
      WomArrayType(WomStringType)
  }

  it should "convert Cwl Input String" in {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    Coproduct[MyriadInputType](y).fold(MyriadInputTypeToWomType).apply(SchemaDefRequirement()) shouldBe WomStringType
  }

  it should "convert Array of a single type is actually one type not in an array" in {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    Coproduct[MyriadInputType](Array(y)).fold(MyriadInputTypeToWomType).apply(SchemaDefRequirement()) shouldBe
      WomStringType
  }

  it should "convert Array of more than one type becomes a coproduct" in {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    val z = Coproduct[MyriadInputInnerType](CwlType.Boolean)
    Coproduct[MyriadInputType](Array(y, z)).fold(MyriadInputTypeToWomType).apply(SchemaDefRequirement()) shouldBe
      WomCoproductType(NonEmptyList.of(WomStringType, WomBooleanType))
  }

  it should "convert Array of more than one type and a null becomes an optional coproduct" in {
    val x = Coproduct[MyriadInputInnerType](CwlType.Null)
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    val z = Coproduct[MyriadInputInnerType](CwlType.Boolean)
    Coproduct[MyriadInputType](Array(x, y, z)).fold(MyriadInputTypeToWomType).apply(SchemaDefRequirement()) shouldBe
      WomOptionalType(WomCoproductType(NonEmptyList.of(WomStringType, WomBooleanType)))
  }

  it should "convert a 2-element Input Array of a single type accompanied by a null is an optional type" in {
    val y = Coproduct[MyriadInputInnerType](CwlType.String)
    val z = Coproduct[MyriadInputInnerType](CwlType.Null)
    testInputArray(Array(y,z), WomOptionalType(WomStringType))
  }

  private def testInputArray(array: Array[MyriadInputInnerType], assertedType: WomType ) = {
    def f(array: Array[MyriadInputInnerType]) =
      Try(Coproduct[MyriadInputType](array).fold(MyriadInputTypeToWomType).apply(SchemaDefRequirement())).cata(
        success => withClue("input should evaluate to an optional array of files type")(success shouldBe assertedType),
        failure => fail(s"expected an optional array of files type but received a failure: $failure")
      )
    f(array)
    f(array.reverse)
  }

  it should "convert Optional Array of a type is interpreted correctly as input type" in {
    val miit = Coproduct[MyriadInputInnerType](CwlType.File)
    val mit = Coproduct[MyriadInputType](miit)
    val ias =  InputArraySchema(items = mit)
    val y = Coproduct[MyriadInputInnerType](ias)
    val z = Coproduct[MyriadInputInnerType](CwlType.Null)
    testInputArray(Array(y, z), WomOptionalType(WomMaybeEmptyArrayType(WomMaybePopulatedFileType)))
  }

  /* ******** Outputs *********** */
  it should "convert ArrayOutputSchema" in {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    val x = OutputArraySchema(items = Coproduct[MyriadOutputType](y))
    val z = Coproduct[MyriadOutputInnerType](x)
    Coproduct[MyriadOutputType](z).fold(MyriadOutputTypeToWomType).apply(SchemaDefRequirement()) shouldBe
      WomArrayType(WomStringType)
  }

  it should "convert Cwl Output String" in {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    Coproduct[MyriadOutputType](y).fold(MyriadOutputTypeToWomType).apply(SchemaDefRequirement()) shouldBe WomStringType
  }

  it should "convert Array of a single type is the same as one type not in an array" in {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    Coproduct[MyriadOutputType](Array(y)).fold(MyriadOutputTypeToWomType).apply(SchemaDefRequirement()) shouldBe
      WomStringType
  }

  it should "convert Output Array of more than one type becomes a coproduct" in {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    val z = Coproduct[MyriadOutputInnerType](CwlType.Boolean)
    Coproduct[MyriadOutputType](Array(y, z)).fold(MyriadOutputTypeToWomType).apply(SchemaDefRequirement()) shouldBe
      WomCoproductType(NonEmptyList.of(WomStringType, WomBooleanType))
  }

  it should "convert Output Array of more than one types including a null is an optional coproduct" in {
    val x = Coproduct[MyriadOutputInnerType](CwlType.Null)
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    val z = Coproduct[MyriadOutputInnerType](CwlType.Boolean)
    Coproduct[MyriadOutputType](Array(x, y, z)).fold(MyriadOutputTypeToWomType).apply(SchemaDefRequirement()) shouldBe
      WomOptionalType(WomCoproductType(NonEmptyList.of(WomStringType, WomBooleanType)))
  }

  private def testOutputArray(array: Array[MyriadOutputInnerType], assertedType: WomType) = {
    def f(array: Array[MyriadOutputInnerType]) =
      Try(Coproduct[MyriadOutputType](array).fold(MyriadOutputTypeToWomType).apply(SchemaDefRequirement())).cata(
        success => withClue("input should evaluate to an optional string type")(success shouldBe assertedType),
        failure => fail(s"expected an optional string type but received a failure: $failure")
      )
    f(array)
    f(array.reverse)
  }

  it should "convert a 2-element Output Array of a single type accompanied by a null is an optional type" in {
    val y = Coproduct[MyriadOutputInnerType](CwlType.String)
    val z = Coproduct[MyriadOutputInnerType](CwlType.Null)
    testOutputArray(Array(y,z), WomOptionalType(WomStringType))
  }

  it should "convert Optional Array of a type is interpreted correctly as output type" in {
    val moit = Coproduct[MyriadOutputInnerType](CwlType.File)
    val mot = Coproduct[MyriadOutputType](moit)
    val oas =  OutputArraySchema(items = mot)
    val moit2 = Coproduct[MyriadOutputInnerType](oas)
    val moit3 = Coproduct[MyriadOutputInnerType](CwlType.Null)
    testOutputArray(Array(moit2,moit3), WomOptionalType(WomMaybeEmptyArrayType(WomMaybePopulatedFileType)))
  }
}
