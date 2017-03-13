package cromwell.backend.standard

import com.typesafe.config.Config
import cromwell.backend.validation._

/**
  * Validates a collection of runtime attributes for a Standard Backend.
  *
  * There are two collections of runtime attributes, the first that are required for the standard async execution, and
  * custom validations that the backend sub class may specify.
  *
  * Currently the required validations are:
  * - `ContinueOnReturnCodeValidation.default`
  * - `FailOnStderrValidation.default`
  *
  * NOTE: The required runtime attributes may be moved to the engine in the future.
  */
object StandardValidatedRuntimeAttributesBuilder {

  private case class StandardValidatedRuntimeAttributesBuilderImpl
  (
    override val requiredValidations: Seq[RuntimeAttributesValidation[_]],
    override val customValidations: Seq[RuntimeAttributesValidation[_]]
  ) extends StandardValidatedRuntimeAttributesBuilder

  /**
    * `default` returns the default set of attributes required to run a standard backend:
    * - `ContinueOnReturnCodeValidation.default`
    * - `FailOnStderrValidation.default`
    *
    * Additional runtime attribute validations may be added by calling `withValidation` on the default.
    */
  def default(backendRuntimeConfig: Config): StandardValidatedRuntimeAttributesBuilder = {
    val required = Seq(ContinueOnReturnCodeValidation.default(backendRuntimeConfig), FailOnStderrValidation.default(backendRuntimeConfig))
    val custom = Seq.empty
    StandardValidatedRuntimeAttributesBuilderImpl(custom, required)
  }

  private def withValidations(builder: StandardValidatedRuntimeAttributesBuilder,
                              customValidations: Seq[RuntimeAttributesValidation[_]]):
  StandardValidatedRuntimeAttributesBuilder = {
    val required = builder.requiredValidations
    val custom = builder.customValidations ++ customValidations
    StandardValidatedRuntimeAttributesBuilderImpl(custom, required)
  }
}

sealed trait StandardValidatedRuntimeAttributesBuilder extends ValidatedRuntimeAttributesBuilder {
  /**
    * Returns a new builder with the additional validation(s).
    *
    * @param validation Additional validation.
    * @return New builder with the validation.
    */
  final def withValidation(validation: RuntimeAttributesValidation[_]*):
  StandardValidatedRuntimeAttributesBuilder = {
    StandardValidatedRuntimeAttributesBuilder.withValidations(this, validation)
  }

  /** Returns all the validations, those required for the standard backend, plus custom addons for the subclass. */
  override final lazy val validations: Seq[RuntimeAttributesValidation[_]] = requiredValidations ++ customValidations

  private[standard] def requiredValidations: Seq[RuntimeAttributesValidation[_]]

  private[standard] def customValidations: Seq[RuntimeAttributesValidation[_]]
}
