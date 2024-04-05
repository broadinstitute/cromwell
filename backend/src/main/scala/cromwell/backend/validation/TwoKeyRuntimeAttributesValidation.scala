package cromwell.backend.validation

import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.RuntimeAttributeDefinition
import wom.values.{WomString, WomValue}

trait TwoKeyRuntimeAttributesValidation[A] extends RuntimeAttributesValidation[A] {

  /**
   * Returns the alternate key of the runtime attribute.
   *
   * @return the alternate key of the runtime attribute.
   */
  def altKey: String

  /**
   * The default value of this runtime attribute.
   *
   * @return the default value of this runtime attribute.
   */
  def defaultVal: A

  /**
   * Runs this validation on the value matching key.
   *
   * NOTE: The values passed to this method should already be evaluated instances of WomValue, and not WomExpression.
   *
   * @param values The full set of values.
   * @return The error or valid value for this key.
   */
  def validateAltKey(values: Map[String, WomValue]): ErrorOr[A] =
    values.get(altKey) match {
      case Some(value) => validateValue.applyOrElse(value, (_: Any) => invalidValueFailure(value))
      case None => validateNone
    }

  /**
   * Returns the value of the default runtime attribute of a
   * validation key as specified in the reference.conf. Given
   * a value, this method coerces it into an optional
   * WomValue. In case the value cannot be succesfully coerced
   * the value is wrapped as a "BadDefaultAttributeValue" type that
   * is failed downstream by the ValidatedRuntimeAttributesBuilder.
   *
   * @param optionalRuntimeConfig Optional default runtime attributes config of a particular backend.
   * @return The new version of this validation.
   */
  final def configDefault(optionalRuntimeConfig: Option[Config]): Option[WomValue] =
    optionalRuntimeConfig collect {
      case config if config.hasPath(key) || config.hasPath(altKey) =>
        val value = if (config.hasPath(key)) config.getValue(key).unwrapped() else config.getValue(altKey).unwrapped()
        coercion collectFirst {
          case womType if womType.coerceRawValue(value).isSuccess => womType.coerceRawValue(value).get
        } getOrElse {
          BadDefaultAttribute(WomString(value.toString))
        }
    }

  /**
   * Returns this as an instance of a runtime attribute definition.
   */
  final lazy val altKeyRuntimeAttributeDefinition =
    RuntimeAttributeDefinition(altKey, staticDefaultOption, usedInCallCaching)
}
