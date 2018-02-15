package wom.types

import org.scalatest.{FlatSpec, Matchers}
import wom.values.{WomBoolean, WomCoproductValue}

class WomCoproductSpec extends FlatSpec with Matchers {

  behavior of "Wom Coproduct Type"

  it should "pick the exact type to coerce to" in {
    val wct = WomCoproductType(List(WomStringType, WomBooleanType))

    wct.coercion(WomBoolean(true)) shouldBe WomCoproductValue(wct, WomBoolean(true))
  }

}
