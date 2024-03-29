package cromwell.backend.validation

import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.RuntimeAttributeDefinition
import cromwell.backend.validation.RuntimeAttributesValidation.extractOption
import wom.types.WomType
import wom.values.{WomString, WomValue}

object TwoKeyRuntimeAttributesValidation {
  def withDefault[ValidatedType, DefaultValType](
    validation: TwoKeyRuntimeAttributesValidation[ValidatedType, DefaultValType],
    default: WomValue
  ): TwoKeyRuntimeAttributesValidation[ValidatedType, DefaultValType] =
    new TwoKeyRuntimeAttributesValidation[ValidatedType, DefaultValType] {
      override def key: String = validation.key

      override def altKey: String = validation.altKey

      override def defaultVal: DefaultValType = validation.defaultVal

      override def coercion: Iterable[WomType] = validation.coercion

      override protected def validateValue: PartialFunction[WomValue, ErrorOr[ValidatedType]] =
        validation.validateValuePackagePrivate

      override protected def validateExpression: PartialFunction[WomValue, Boolean] =
        validation.validateExpressionPackagePrivate

      override protected def invalidValueMessage(value: WomValue): String =
        validation.invalidValueMessagePackagePrivate(value)

      override protected def missingValueMessage: String = s"Expecting $key runtime attribute to be a type in $coercion"

      override def usedInCallCaching: Boolean = validation.usedInCallCachingPackagePrivate

      override protected def staticDefaultOption = Option(default)
    }

  /**
   * Returns the value from the attributes matching one of the validation keys. If values are assigned to attribute
   * matching both keys, the `key` field will override the `altKey` field on the `RuntimeAttributeValidation`.
   *
   * Do not use an optional validation as the type internal implementation will throw a `ClassCastException` due to the
   * way values are located and auto-magically cast to the type of the `runtimeAttributesValidation`.
   *
   * @param runtimeAttributesValidation The typed validation to use.
   * @param validatedRuntimeAttributes  The values to search.
   * @return The value matching the key.
   * @throws ClassCastException if the validation is called on an optional validation.
   */
  def extractTwoKeys[A, B](runtimeAttributesValidation: TwoKeyRuntimeAttributesValidation[A, B],
                           validatedRuntimeAttributes: ValidatedRuntimeAttributes
  ): A = {
    val key = runtimeAttributesValidation.key
    val altKey = runtimeAttributesValidation.altKey
    val primaryValue = extractOption(key, validatedRuntimeAttributes)
    val altValue = extractOption(altKey, validatedRuntimeAttributes)
    val value =
      if (
        primaryValue.isEmpty || (primaryValue.contains(
          runtimeAttributesValidation.defaultVal
        ) && !altValue.isEmpty && !altValue.contains(
          runtimeAttributesValidation.defaultVal
        ))
      ) altValue
      else primaryValue
    value match {
      // NOTE: Some(innerValue) aka Some.unapply() throws a `ClassCastException` to `Nothing$` as it can't tell the type
      case some: Some[_] => some.get.asInstanceOf[A]
      case None =>
        throw new RuntimeException(
          s"$key not found in runtime attributes ${validatedRuntimeAttributes.attributes.keys}"
        )
    }
  }
}

trait TwoKeyRuntimeAttributesValidation[A, DefaultValType] extends RuntimeAttributesValidation[A] {

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
  def defaultVal: DefaultValType

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
   * Returns a version of this validation with the default value.
   *
   * @param womValue The default wdl value.
   * @return The new version of this validation.
   */
  final def makeDefault(womValue: WomValue): TwoKeyRuntimeAttributesValidation[A, DefaultValType] =
    TwoKeyRuntimeAttributesValidation.withDefault(this, womValue)

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
