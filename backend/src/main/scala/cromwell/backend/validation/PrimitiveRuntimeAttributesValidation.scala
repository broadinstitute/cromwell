package cromwell.backend.validation

import cats.syntax.validated._
import lenthall.validation.ErrorOr.ErrorOr
import wom.types._
import wom.values._

/**
  * Validates one of the wdl primitive types: Boolean, Float, Integer, or String. WdlFile is not supported.
  *
  * @tparam A The type of validated runtime attribute.
  * @tparam B The type of coerced WomValue.
  */
sealed trait PrimitiveRuntimeAttributesValidation[A, B <: WomPrimitive] extends RuntimeAttributesValidation[A] {
  val womType: WomPrimitiveType

  override def coercion = Seq(womType)

  override protected def validateExpression: PartialFunction[WomValue, Boolean] = {
    case womValue if womType.coerceRawValue(womValue).isSuccess => true
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be $typeString"

  protected def typeString = s"a ${womType.toDisplayString}"

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[A]] = {
    case value if womType.coerceRawValue(value).isSuccess =>
      validateCoercedValue(womType.coerceRawValue(value).get.asInstanceOf[B])
    case value if womType.coerceRawValue(value.valueString).isSuccess =>
      /*
      NOTE: This case statement handles WdlString("true") coercing to WdlBoolean(true).
      For some reason "true" as String is coercable... but not the WdlString.
       */
      validateCoercedValue(womType.coerceRawValue(value.valueString).get.asInstanceOf[B])
  }

  protected def validateCoercedValue(womValue: B): ErrorOr[A]
}

class BooleanRuntimeAttributesValidation(override val key: String) extends
  PrimitiveRuntimeAttributesValidation[Boolean, WomBoolean] {

  override val womType = WomBooleanType

  override protected def validateCoercedValue(womValue: WomBoolean): ErrorOr[Boolean] = womValue.value.validNel
}

class FloatRuntimeAttributesValidation(override val key: String) extends
  PrimitiveRuntimeAttributesValidation[Double, WomFloat] {

  override val womType = WomFloatType

  override protected def validateCoercedValue(womValue: WomFloat): ErrorOr[Double] = womValue.value.validNel
}

class IntRuntimeAttributesValidation(override val key: String) extends
  PrimitiveRuntimeAttributesValidation[Int, WomInteger] {

  override val womType = WomIntegerType

  override protected def validateCoercedValue(womValue: WomInteger): ErrorOr[Int] = womValue.value.toInt.validNel

  override protected def typeString: String = "an Integer"
}

class StringRuntimeAttributesValidation(override val key: String) extends
  PrimitiveRuntimeAttributesValidation[String, WomString] {

  override val womType = WomStringType

  override protected def validateCoercedValue(womValue: WomString): ErrorOr[String] = womValue.value.validNel
}
