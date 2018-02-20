package cwl

import cats.data.NonEmptyList
import org.scalatest.{FlatSpec, Matchers}
import shapeless.Coproduct
import wom.types.{WomBooleanType, WomCoproductType, WomOptionalType, WomStringType}

class MyriadInputToWomTypeSpec extends FlatSpec with Matchers {
  behavior of "Myriad Input To Wom Type"

  it should "convert multiple types into a coproduct type" in {
    val bool = Coproduct[MyriadInputInnerType](CwlType.Boolean)
    val string = Coproduct[MyriadInputInnerType](CwlType.String)
    val array = Array(bool, string)
    val mit = Coproduct[MyriadInputType](array)

    mit.fold(MyriadInputTypeToWomType) shouldBe WomCoproductType(NonEmptyList.of(WomBooleanType, WomStringType))
  }

  it should "convert optional multiple types into an optional coproduct type" in {
    val Null = Coproduct[MyriadInputInnerType](CwlType.Null)
    val bool = Coproduct[MyriadInputInnerType](CwlType.Boolean)
    val string = Coproduct[MyriadInputInnerType](CwlType.String)
    val array = Array(bool, string, Null)
    val mit = Coproduct[MyriadInputType](array)

    mit.fold(MyriadInputTypeToWomType) shouldBe WomOptionalType(WomCoproductType(NonEmptyList.of(WomBooleanType, WomStringType)))
  }

}


class MyriadOutput extends FlatSpec with Matchers {
  behavior of "Myriad Output To Wom Type"

  it should "convert multiple types into a coproduct type" in {
    val bool = Coproduct[MyriadOutputInnerType](CwlType.Boolean)
    val string = Coproduct[MyriadOutputInnerType](CwlType.String)
    val array = Array(bool, string)
    val mit = Coproduct[MyriadOutputType](array)

    mit.fold(MyriadOutputTypeToWomType) shouldBe WomCoproductType(NonEmptyList.of(WomBooleanType, WomStringType))
  }

  it should "convert optional multiple types into an optional coproduct type" in {
    val Null = Coproduct[MyriadOutputInnerType](CwlType.Null)
    val bool = Coproduct[MyriadOutputInnerType](CwlType.Boolean)
    val string = Coproduct[MyriadOutputInnerType](CwlType.String)
    val array = Array(bool, string, Null)
    val mit = Coproduct[MyriadOutputType](array)

    mit.fold(MyriadOutputTypeToWomType) shouldBe WomOptionalType(WomCoproductType(NonEmptyList.of(WomBooleanType, WomStringType)))
  }

}

