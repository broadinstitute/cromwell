package wom.types

import wom.values.{WomInteger, WomLong, WomString}
import org.scalacheck.Properties
import org.scalacheck.Prop._
import spray.json.{JsNumber, JsString}

import scala.util.Success


object WomLongTypeSpec extends Properties("WomLongType") {

  property("conversion from Long") = forAll { i: Long =>
    WomLongType.coerceRawValue(i) == Success(WomLong(i))
  }

  property("conversion from String") = forAll { i: Long =>
    WomLongType.coerceRawValue(i.toString) == Success(WomLong(i))
  }

  property("conversion from Wom String") = forAll { i: Long =>
    WomLongType.coerceRawValue(WomString(i.toString)) == Success(WomLong(i))
  }

  property("conversion from Js String") = forAll { i: Long =>
    WomLongType.coerceRawValue(JsString(i.toString)) == Success(WomLong(i))
  }

  property("conversion from Int") = forAll { i: Int =>
    WomLongType.coerceRawValue(i) == Success(WomLong(i.toLong))
  }

  property("conversion from WomInt") = forAll { i: Int =>
    WomLongType.coerceRawValue(WomInteger(i)) == Success(WomLong(i.toLong))
  }

  property("conversion from JsNumber") = forAll { i: Long =>
    WomLongType.coerceRawValue(JsNumber(i)) == Success(WomLong(i))
  }
}

