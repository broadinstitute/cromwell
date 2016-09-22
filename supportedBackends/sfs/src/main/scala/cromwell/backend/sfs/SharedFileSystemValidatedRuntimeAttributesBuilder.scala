package cromwell.backend.sfs

import cromwell.backend.validation._

/**
  * Validates a collection of runtime attributes for a Shared File System (SFS) Backend.
  *
  * There are always three collections of runtime attributes, two of which are required for the SFS base class:
  *
  * Required for the SFS to be able to run properly:
  * 1) The set of runtime attributes that are absolutely required, or the validation fails.
  * 2) _Extra_ validations that will run for the SFS.
  *
  * Lastly for the sub classes of the SFS, there are:
  * 3) Custom validations that the backend sub class may specify.
  *
  * 3), the custom validations, are always set by the sub class, via calls to `withValidation()`.
  *
  * For 1) and 2) the biggest difference is when it comes to docker support in a backend.
  *
  * For a backend that does __not__ support docker, the default, 1) and 2) will contain:
  *
  * 1) required         = ContinueOnReturnCodeValidation.default, FailOnStderrValidation.default
  * 2) unsupportedExtra = DockerValidation.optional
  *
  * Suppose the above validation runs a WDL with runtime attributes shows up with:
  *
  * {{{
  *   runtimeAttributes {
  *     continueOnReturnCode: 14
  *     # failOnStdErr not specified
  *     docker: "ubuntu"
  *   }
  * }}}
  *
  * This will cause a validation warning to print out that docker is unsupported, because the docker validation isn't
  * present in the required validations. If/when the SFS backend asks what the value of the optional docker is, it will
  * receive `None`, as the validation _is not_ listed in the required validations.
  *
  * There is an additional interesting thing about having the docker validation still running under an "unsupported
  * extra validation". Say an invalid docker were to be specified as a WdlArray. The extra validation __would__ catch
  * the error, even after the warning had been printed stating that docker is not supported by the backend.
  *
  * Meanwhile, even though failOnStdErr is not specified, the `FailOnStderrValidation.default` will return its default
  * value. And of course, the `ContinueOnReturnCodeValidation.default` returns the specified, and valid, runtime
  * attribute `ContinueOnReturnCodeSet(14)`.
  *
  * Now--
  *
  * Suppose the `withDockerSupport(true)` has been invoked on the builder. The required and unsupported runtime
  * attributes will then look like:
  *
  * 1) required         = ContinueOnReturnCodeValidation.default, FailOnStderrValidation.default,
  *                       DockerValidation.optional
  * 2) unsupportedExtra = __empty__
  *
  * With the same WDL above, this builder does NOT print a warning, because docker __is__ supported. When the SFS asks
  * for the optional docker element, it receives `Some("ubuntu")`, as the `DockerValidation.optional` __is__ listed in
  * the required validations.
  *
  * `ContinueOnReturnCodeValidation.default` and `FailOnStderrValidation.default` still operate as in the previous
  * example.
  *
  * What happens when there is no runtime attribute for docker? Easy, the validation for docker is always optional!
  * In either case of running a builder via `withDockerSupport(true)` or `withDockerSupport(false)`, if the docker
  * runtime attribute is not specified, the `SharedFileSystemValidatedRuntimeAttributesBuilder` will return a
  * `None` value.
  */
object SharedFileSystemValidatedRuntimeAttributesBuilder {

  private case class SharedFileSystemValidatedRuntimeAttributesBuilderImpl
  (
    override val requiredValidations: Seq[RuntimeAttributesValidation[_]],
    override val customValidations: Seq[RuntimeAttributesValidation[_]]
  ) extends SharedFileSystemValidatedRuntimeAttributesBuilder

  /**
    * `default` returns the default set of attributes required to run an SFS:
    * - `ContinueOnReturnCodeValidation.default`
    * - `FailOnStderrValidation.default`
    *
    * Additional runtime attribute validations may be added by calling `withValidation` on the default.
    *
    * The SFS will also _always_ validate using the `DockerValidation`, but will end up warning the user that the
    * runtime attribute is unsupported by the backend implementation unless `withDockerSupport(true)` is called.
    */
  lazy val default: SharedFileSystemValidatedRuntimeAttributesBuilder = {
    val required = Seq(ContinueOnReturnCodeValidation.default, FailOnStderrValidation.default)
    val custom = Seq.empty
    SharedFileSystemValidatedRuntimeAttributesBuilderImpl(custom, required)
  }

  private def withValidations(builder: SharedFileSystemValidatedRuntimeAttributesBuilder,
                              customValidations: Seq[RuntimeAttributesValidation[_]]):
  SharedFileSystemValidatedRuntimeAttributesBuilder = {
    val required = builder.requiredValidations
    val custom = builder.customValidations ++ customValidations
    SharedFileSystemValidatedRuntimeAttributesBuilderImpl(custom, required)
  }
}

sealed trait SharedFileSystemValidatedRuntimeAttributesBuilder extends ValidatedRuntimeAttributesBuilder {
  /**
    * Returns a new builder with the additional validation(s).
    *
    * @param validation Additional validation.
    * @return New builder with the validation.
    */
  final def withValidation(validation: RuntimeAttributesValidation[_]*):
  SharedFileSystemValidatedRuntimeAttributesBuilder = {
    SharedFileSystemValidatedRuntimeAttributesBuilder.withValidations(this, validation)
  }

  /** Returns on the supported validations, those required for the SFS, plus custom addons for the subclass. */
  override final lazy val validations: Seq[RuntimeAttributesValidation[_]] = requiredValidations ++ customValidations

  private[sfs] def requiredValidations: Seq[RuntimeAttributesValidation[_]]

  private[sfs] def customValidations: Seq[RuntimeAttributesValidation[_]]
}
