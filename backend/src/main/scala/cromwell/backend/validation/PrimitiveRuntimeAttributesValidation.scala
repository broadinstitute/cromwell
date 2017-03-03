package cromwell.backend.validation

import cats.syntax.validated._
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.types._
import wdl4s.values.{WdlBoolean, WdlFloat, WdlInteger, WdlPrimitive, WdlString, WdlValue}

/**
  * Validates one of the wdl primitive types: Boolean, Float, Integer, or String. WdlFile is not supported.
  *
  * @tparam A The type of validated runtime attribute.
  * @tparam B The type of coerced WdlValue.
  */
sealed trait PrimitiveRuntimeAttributesValidation[A, B <: WdlPrimitive] extends RuntimeAttributesValidation[A] {
  val wdlType: WdlPrimitiveType

  override def coercion = Seq(wdlType)

  override protected def validateExpression: PartialFunction[WdlValue, Boolean] = {
    case wdlValue if wdlType.coerceRawValue(wdlValue).isSuccess => true
  }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be $typeString"

  protected def typeString = s"a ${wdlType.toWdlString}"

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[A]] = {
    case value if wdlType.coerceRawValue(value).isSuccess =>
      validateCoercedValue(wdlType.coerceRawValue(value).get.asInstanceOf[B])
    case value if wdlType.coerceRawValue(value.valueString).isSuccess =>
      /*
      NOTE: This case statement handles WdlString("true") coercing to WdlBoolean(true).
      For some reason "true" as String is coercable... but not the WdlString.
       */
      validateCoercedValue(wdlType.coerceRawValue(value.valueString).get.asInstanceOf[B])
  }

  protected def validateCoercedValue(wdlValue: B): ErrorOr[A]
}

class BooleanRuntimeAttributesValidation(override val key: String) extends
  PrimitiveRuntimeAttributesValidation[Boolean, WdlBoolean] {

  override val wdlType = WdlBooleanType

  override protected def validateCoercedValue(wdlValue: WdlBoolean): ErrorOr[Boolean] = wdlValue.value.validNel
}

class FloatRuntimeAttributesValidation(override val key: String) extends
  PrimitiveRuntimeAttributesValidation[Double, WdlFloat] {

  override val wdlType = WdlFloatType

  override protected def validateCoercedValue(wdlValue: WdlFloat): ErrorOr[Double] = wdlValue.value.validNel
}

class IntRuntimeAttributesValidation(override val key: String) extends
  PrimitiveRuntimeAttributesValidation[Int, WdlInteger] {

  override val wdlType = WdlIntegerType

  override protected def validateCoercedValue(wdlValue: WdlInteger): ErrorOr[Int] = wdlValue.value.toInt.validNel

  override protected def typeString: String = "an Integer"
}

class StringRuntimeAttributesValidation(override val key: String) extends
  PrimitiveRuntimeAttributesValidation[String, WdlString] {

  override val wdlType = WdlStringType

  override protected def validateCoercedValue(wdlValue: WdlString): ErrorOr[String] = wdlValue.value.validNel
}
