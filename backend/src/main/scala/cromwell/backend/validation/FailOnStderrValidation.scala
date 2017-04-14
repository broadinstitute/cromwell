package cromwell.backend.validation

import com.typesafe.config.Config
import wdl4s.values.{WdlBoolean, WdlValue}

/**
  * Validates the "failOnStderr" runtime attribute as a Boolean or a String 'true' or 'false', returning the value as a
  * `Boolean`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * `configDefaultWdlValue` returns the value of the attribute as specified by the
  * reference.conf file, coerced into a WdlValue.
  *
  * `default` a validation with the default value specified by the reference.conf file.
  */

object FailOnStderrValidation {
  lazy val instance: RuntimeAttributesValidation[Boolean] = new FailOnStderrValidation
  def default(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Boolean] = instance.withDefault(
    configDefaultWdlValue(runtimeConfig) getOrElse WdlBoolean(false))
  def configDefaultWdlValue(runtimeConfig: Option[Config]): Option[WdlValue] = instance.configDefaultWdlValue(runtimeConfig)
}

class FailOnStderrValidation extends BooleanRuntimeAttributesValidation(RuntimeAttributesKeys.FailOnStderrKey) {
  override protected def usedInCallCaching: Boolean = true

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be a Boolean or a String with values of 'true' or 'false'"
}
