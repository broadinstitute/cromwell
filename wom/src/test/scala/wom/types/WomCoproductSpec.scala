package wom.types

import cats.data.NonEmptyList
import org.scalatest.{FlatSpec, Matchers}
import wom.values.{WomBoolean, WomCoproductValue, WomString}

class WomCoproductSpec extends FlatSpec with Matchers {

  behavior of "Wom Coproduct Type"

  it should "pick the exact type to coerce to" in {
    val wct = WomCoproductType(NonEmptyList.of(WomStringType, WomBooleanType))

    wct.coercion(WomBoolean(true)) shouldBe WomCoproductValue(wct, WomBoolean(true))
  }

  it should "attempt to coerce to expected types in order" in {
    val wct = WomCoproductType(NonEmptyList.of(WomStringType, WomBooleanType))

    wct.coercion(true) shouldBe WomCoproductValue(wct, WomBoolean(true))
    wct.coercion("true") shouldBe WomCoproductValue(wct, WomString("true"))
  }
}
