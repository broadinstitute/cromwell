package wom.types

import org.scalatest.prop.{TableDrivenPropertyChecks, TableFor2}
import org.scalatest.{FlatSpecLike, Matchers}
import wom.values.WomValue

import scala.util.{Failure, Success}

abstract class WomCoercionSpec(val goodCoercionTable: TableFor2[_ <: Any, WomValue],
                               val badCoercionTable: TableFor2[_ <: Any, WomType],
                               val behaviorOf: String) extends FlatSpecLike with Matchers {

  import TableDrivenPropertyChecks._
  import WomCoercionSpec.StringableAny

  behavior of behaviorOf

  forAll(goodCoercionTable) { (fromValue, toValue) =>

    it should s"Allow coercion from ${fromValue.displayString} (${fromValue.typeDisplayString}) to ${toValue.toWomString} (${toValue.womType.stableName})" in {
      toValue.womType.coercionDefined(fromValue) should be(true)
      fromValue match {
        case wv: WomValue => 
          toValue.womType.isCoerceableFrom(wv.womType) should be(true)
        case _ => // can't test isCoerceableFrom for this fromValue
      }
    }

    it should s"correctly coerce ${fromValue.displayString} (${fromValue.typeDisplayString}) to ${toValue.toWomString} (${toValue.womType.stableName})" in {
      toValue.womType.coerceRawValue(fromValue) should be(Success(toValue))
    }
  }

  forAll(badCoercionTable) { (fromValue, toType) =>

    it should s"Not allow coercion from ${fromValue.displayString} to ${toType.stableName}" in {
      toType.coercionDefined(fromValue) should be(false)
      fromValue match {
        case wv: WomValue => toType.isCoerceableFrom(wv.womType) should be(false)
        case _ => // can't test isCoerceableFrom for this fromValue
      }
    }

    it should s"Fail when coercing ${fromValue.displayString} to ${toType.stableName}" in {
      toType.coerceRawValue(fromValue) match {
        case Failure(_) => // great!
        case Success(v) => fail(s"Incorrectly coerced $fromValue to $v")
      }
    }
  }
}

private object WomCoercionSpec {
  implicit class StringableAny(val a: Any) extends AnyVal {
    def displayString: String = a match {
      case wv: WomValue => wv.toWomString
      case other => other.toString
    }

    def typeDisplayString: String = a match {
      case wv: WomValue => wv.womType.stableName
      case other => other.getClass.getSimpleName
    }
  }
}
