package cromwell.backend.impl.sfs.config

import cats.syntax.validated._
import cromwell.backend.validation.RuntimeAttributesValidation
import wdl4s.types._
import wdl4s.values.{WdlBoolean, WdlFloat, WdlInteger, WdlString}

/**
  * Validates one of the wdl primitive types: Boolean, Float, Integer, or String. WdlFile is not supported.
  *
  * @tparam A The type of validated runtime attribute.
  */
sealed trait PrimitiveRuntimeAttributesValidation[A] extends RuntimeAttributesValidation[A] {
  val wdlType: WdlPrimitiveType

  override def coercion = Seq(wdlType)
}

class BooleanRuntimeAttributesValidation(override val key: String) extends
  PrimitiveRuntimeAttributesValidation[Boolean] {

  override val wdlType = WdlBooleanType

  override protected def validateValue = {
    case WdlBoolean(value) => value.validNel
  }
}

class FloatRuntimeAttributesValidation(override val key: String) extends PrimitiveRuntimeAttributesValidation[Double] {
  override val wdlType = WdlFloatType

  override protected def validateValue = {
    case WdlFloat(value) => value.validNel
  }
}

class IntRuntimeAttributesValidation(override val key: String) extends PrimitiveRuntimeAttributesValidation[Int] {
  override val wdlType = WdlIntegerType

  override protected def validateValue = {
    case WdlInteger(value) => value.toInt.validNel
  }
}

class StringRuntimeAttributesValidation(override val key: String) extends PrimitiveRuntimeAttributesValidation[String] {
  override val wdlType = WdlStringType

  override protected def validateValue = {
    case WdlString(value) => value.validNel
  }
}
