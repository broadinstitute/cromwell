package wom.types

import wom.values.{WomInteger, WomLong, WomString}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatestplus.scalacheck.ScalaCheckDrivenPropertyChecks
import spray.json.{JsNumber, JsString}

import scala.util.Success

class WomLongTypeSpec extends AnyFlatSpec with Matchers with ScalaCheckDrivenPropertyChecks {
  behavior of "WomLongType"

  it should "conversion from Long" in forAll { i: Long =>
    WomLongType.coerceRawValue(i) shouldBe Success(WomLong(i))
  }

  it should "conversion from String" in forAll { i: Long =>
    WomLongType.coerceRawValue(i.toString) shouldBe Success(WomLong(i))
  }

  it should "conversion from Wom String" in forAll { i: Long =>
    WomLongType.coerceRawValue(WomString(i.toString)) shouldBe Success(WomLong(i))
  }

  it should "conversion from Js String" in forAll { i: Long =>
    WomLongType.coerceRawValue(JsString(i.toString)) shouldBe Success(WomLong(i))
  }

  it should "conversion from Int" in forAll { i: Int =>
    WomLongType.coerceRawValue(i) shouldBe Success(WomLong(i.toLong))
  }

  it should "conversion from WomInt" in forAll { i: Int =>
    WomLongType.coerceRawValue(WomInteger(i)) shouldBe Success(WomLong(i.toLong))
  }

  it should "conversion from JsNumber" in forAll { i: Long =>
    WomLongType.coerceRawValue(JsNumber(i)) shouldBe Success(WomLong(i))
  }
}
