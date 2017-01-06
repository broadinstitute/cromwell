package cromwell.backend.validation

import wdl4s.values.WdlBoolean

/**
  * Validates the "failOnStderr" runtime attribute as a Boolean or a String 'true' or 'false', returning the value as a
  * `Boolean`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * The default returns `false` when no attribute is specified.
  */
object FailOnStderrValidation {
  lazy val instance: RuntimeAttributesValidation[Boolean] = new FailOnStderrValidation
  lazy val default: RuntimeAttributesValidation[Boolean] = instance.withDefault(WdlBoolean(false))
  lazy val optional: OptionalRuntimeAttributesValidation[Boolean] = default.optional
}

class FailOnStderrValidation extends BooleanRuntimeAttributesValidation(RuntimeAttributesKeys.FailOnStderrKey) {
  override protected def usedInCallCaching: Boolean = true

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be a Boolean or a String with values of 'true' or 'false'"
}
