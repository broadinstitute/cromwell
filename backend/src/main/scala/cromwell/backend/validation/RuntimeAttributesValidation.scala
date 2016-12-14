package cromwell.backend.validation

import cats.syntax.validated._
import wdl4s.expression.PureStandardLibraryFunctions
import cromwell.backend.{MemorySize, RuntimeAttributeDefinition}
import lenthall.validation.ErrorOr._
import org.slf4j.Logger
import wdl4s.WdlExpression
import wdl4s.WdlExpression._
import wdl4s.types.{WdlBooleanType, WdlIntegerType, WdlType}
import wdl4s.values._

import scala.util.{Failure, Success}

object RuntimeAttributesValidation {

  def warnUnrecognized(actual: Set[String], expected: Set[String], logger: Logger) = {
    val unrecognized = actual.diff(expected).mkString(", ")
    if (unrecognized.nonEmpty) logger.warn(s"Unrecognized runtime attribute keys: $unrecognized")
  }

  def validateDocker(docker: Option[WdlValue], onMissingKey: => ErrorOr[Option[String]]): ErrorOr[Option[String]] = {
    validateWithValidation(docker, DockerValidation.optional, onMissingKey, DockerValidation.missingMessage)
  }

  def validateFailOnStderr(value: Option[WdlValue], onMissingKey: => ErrorOr[Boolean]): ErrorOr[Boolean] = {
    validateWithValidation(value, FailOnStderrValidation.default, onMissingKey, FailOnStderrValidation.missingMessage)
  }

  def validateContinueOnReturnCode(value: Option[WdlValue],
                                   onMissingKey: => ErrorOr[ContinueOnReturnCode]): ErrorOr[ContinueOnReturnCode] = {
    validateWithValidation(value, ContinueOnReturnCodeValidation.default, onMissingKey,
      ContinueOnReturnCodeValidation.missingMessage)
  }

  def validateMemory(value: Option[WdlValue], onMissingKey: => ErrorOr[MemorySize]): ErrorOr[MemorySize] = {
    validateWithValidation(value, MemoryValidation.instance, onMissingKey,
      MemoryValidation.missingFormat.format("Not supported WDL type value"))
  }

  def validateCpu(cpu: Option[WdlValue], onMissingKey: => ErrorOr[Int]): ErrorOr[Int] = {
    validateWithValidation(cpu, CpuValidation.default, onMissingKey, CpuValidation.missingMessage)
  }

  private def validateWithValidation[T](valueOption: Option[WdlValue],
                                        validation: RuntimeAttributesValidation[T],
                                        onMissingValue: => ErrorOr[T],
                                        missingValidationMessage: String): ErrorOr[T] = {
    valueOption match {
      case Some(value) =>
        validation.validateValue.applyOrElse(value, (_: Any) => missingValidationMessage.invalidNel)
      case None => onMissingValue
    }
  }

  def validateInt(value: WdlValue): ErrorOr[Int] = {
    WdlIntegerType.coerceRawValue(value) match {
      case scala.util.Success(WdlInteger(i)) => i.intValue.validNel
      case _ => s"Could not coerce ${value.valueString} into an integer".invalidNel
    }
  }

  def validateBoolean(value: WdlValue): ErrorOr[Boolean] = {
    WdlBooleanType.coerceRawValue(value) match {
      case scala.util.Success(WdlBoolean(b)) => b.booleanValue.validNel
      case _ => s"Could not coerce ${value.valueString} into a boolean".invalidNel
    }
  }

  def parseMemoryString(s: WdlString): ErrorOr[MemorySize] = {
    MemoryValidation.validateMemoryString(s)
  }

  def parseMemoryInteger(i: WdlInteger): ErrorOr[MemorySize] = {
    MemoryValidation.validateMemoryInteger(i)
  }

  def withDefault[ValidatedType](validation: RuntimeAttributesValidation[ValidatedType],
                                 default: WdlValue): RuntimeAttributesValidation[ValidatedType] = {
    new RuntimeAttributesValidation[ValidatedType] {
      override def key = validation.key

      override def coercion = validation.coercion

      override protected def validateValue = validation.validateValuePackagePrivate

      override protected def validateExpression = validation.validateExpressionPackagePrivate

      override def staticDefaultOption = Option(default)

      override protected def failureMessage = validation.failureMessagePackagePrivate
    }
  }

  def optional[ValidatedType](validation: RuntimeAttributesValidation[ValidatedType]):
  OptionalRuntimeAttributesValidation[ValidatedType] = {
    new OptionalRuntimeAttributesValidation[ValidatedType] {
      override def key = validation.key

      override def coercion = validation.coercion

      override protected def validateOption = validation.validateValuePackagePrivate

      override protected def validateExpression = validation.validateExpressionPackagePrivate

      override protected def failureMessage = validation.failureMessagePackagePrivate
    }
  }

  /**
    * Returns the value from the attributes, unpacking options.
    *
    * @param validatedRuntimeAttributes The values to search.
    * @return The keys and extracted values.
    */
  def extract(validatedRuntimeAttributes: ValidatedRuntimeAttributes): Map[String, Any] = {
    val attributeOptions: Map[String, Option[Any]] = validatedRuntimeAttributes.attributes.mapValues(unpackOption)

    val attributes = attributeOptions collect {
      case (name, Some(value)) => (name, value)
    }

    attributes
  }

  /**
    * Returns the value from the attributes matching the validation key.
    *
    * Do not use an optional validation as the type internal implementation will throw a `ClassCastException` due to the
    * way values are located and auto-magically cast to the type of the `runtimeAttributesValidation`.
    *
    * @param runtimeAttributesValidation The typed validation to use.
    * @param validatedRuntimeAttributes  The values to search.
    * @return The value matching the key.
    * @throws ClassCastException if the validation is called on an optional validation.
    */
  def extract[A](runtimeAttributesValidation: RuntimeAttributesValidation[A],
                 validatedRuntimeAttributes: ValidatedRuntimeAttributes): A = {
    extract(runtimeAttributesValidation.key, validatedRuntimeAttributes)
  }

  /**
    * Returns the value from the attributes matching the key.
    *
    * @param key                        The key to retrieve.
    * @param validatedRuntimeAttributes The values to search.
    * @return The value matching the key.
    */
  def extract[A](key: String,
                 validatedRuntimeAttributes: ValidatedRuntimeAttributes): A = {
    val value = extractOption(key, validatedRuntimeAttributes)
    value match {
      // NOTE: Some(innerValue) aka Some.unapply() throws a `ClassCastException` to `Nothing$` as it can't tell the type
      case some: Some[_] => some.get.asInstanceOf[A]
      case None => throw new RuntimeException(
        s"$key not found in runtime attributes ${validatedRuntimeAttributes.attributes.keys}")
    }
  }

  /**
    * Returns Some(value) from the attributes matching the validation key, or None.
    *
    * @param runtimeAttributesValidation The typed validation to use.
    * @param validatedRuntimeAttributes  The values to search.
    * @return The Some(value) matching the key or None.
    */
  def extractOption[A](runtimeAttributesValidation: RuntimeAttributesValidation[A],
                       validatedRuntimeAttributes: ValidatedRuntimeAttributes): Option[A] = {
    extractOption(runtimeAttributesValidation.key, validatedRuntimeAttributes)
  }

  /**
    * Returns Some(value) from the attributes matching the key, or None.
    *
    * @param key                        The key to retrieve.
    * @param validatedRuntimeAttributes The values to search.
    * @return The Some(value) matching the key or None.
    */
  def extractOption[A](key: String, validatedRuntimeAttributes: ValidatedRuntimeAttributes): Option[A] = {
    val value = validatedRuntimeAttributes.attributes.get(key)
    unpackOption[A](value)
  }

  /**
    * Recursively unpacks an option looking for a value of some type A.
    *
    * @param value The value to unpack.
    * @tparam A The type to cast the unpacked value.
    * @return The Some(value) matching the key or None.
    */
  final def unpackOption[A](value: Any): Option[A] = {
    value match {
      case None => None
      case Some(innerValue) => unpackOption(innerValue)
      case _ => Option(value.asInstanceOf[A])
    }
  }

  /**
    * Converts a RuntimeAttributesValidation to a RuntimeAttributeDefinition.
    *
    * @param validation RuntimeAttributesValidation
    * @return RuntimeAttributeDefinition
    */
  def toRuntimeAttributeDefinition(validation: RuntimeAttributesValidation[_]): RuntimeAttributeDefinition = {
    val name = validation.key
    val default = validation.staticDefaultOption
    import cromwell.backend.validation.RuntimeAttributesKeys._
    val usedInCallCaching = name match {
      case DockerKey | ContinueOnReturnCodeKey | FailOnStderrKey => true
      case _ => false
    }
    RuntimeAttributeDefinition(name, default, usedInCallCaching)
  }
}

/**
  * Performs a validation on a runtime attribute and returns some value.
  *
  * @tparam ValidatedType The type of the validated value.
  */
trait RuntimeAttributesValidation[ValidatedType] {
  /**
    * Returns the key of the runtime attribute.
    *
    * @return The key of the runtime attribute.
    */
  def key: String

  /**
    * The WDL types that will be passed to `validate`, after the value is coerced from the first element found that
    * can coerce the type.
    *
    * @return traversable of wdl types
    */
  def coercion: Traversable[WdlType]

  /**
    * Validates the wdl value.
    *
    * @return The validated value or an error, wrapped in a cats validation.
    */
  protected def validateValue: PartialFunction[WdlValue, ErrorOr[ValidatedType]]

  /**
    * Returns the value for when there is no wdl value. By default returns an error.
    *
    * @return the value for when there is no wdl value.
    */
  protected def validateNone: ErrorOr[ValidatedType] = failureWithMessage

  /**
    * Returns true if the value can be validated.
    *
    * The base implementation does a basic check that a coercion exists.
    *
    * Subclasses may inspect the wdl value for more information to identify if the value may be validated. For example,
    * the `ContinueOnReturnCodeValidation` checks that all elements in a `WdlArray` can be sub-coerced into an integer.
    *
    * @return true if the value can be validated.
    */
  protected def validateExpression: PartialFunction[WdlValue, Boolean] = {
    case wdlValue => coercion.exists(_ == wdlValue.wdlType)
  }

  /**
    * Returns the optional default value when no other is specified.
    *
    * @return the optional default value when no other is specified.
    */
  def staticDefaultOption: Option[WdlValue] = None

  /**
    * Returns message to return when a value is invalid.
    *
    * @return Message to return when a value is invalid.
    */
  protected def failureMessage: String = s"Expecting $key runtime attribute to be a type in $coercion"

  /**
    * Utility method to wrap the failureMessage in an ErrorOr.
    *
    * @return Wrapped failureMessage.
    */
  protected final lazy val failureWithMessage: ErrorOr[ValidatedType] = failureMessage.invalidNel

  /**
    * Runs this validation on the value matching key.
    *
    * @param values The full set of values.
    * @return The error or valid value for this key.
    */
  def validate(values: Map[String, WdlValue]): ErrorOr[ValidatedType] = {
    values.get(key) match {
      case Some(value) => validateValue.applyOrElse(value, (_: Any) => failureWithMessage)
      case None => validateNone
    }
  }

  /**
    * Used during initialization, returning true if the expression __may be__ valid.
    *
    * The `BackendWorkflowInitializationActor` requires validation be performed by a map of validating functions:
    *
    * {{{
    * runtimeAttributeValidators: Map[String, Option[WdlExpression] => Boolean]
    * }}}
    *
    * With our `key` as the key in the map, one can return this function as the value in the map.
    *
    * @param wdlExpressionMaybe The optional expression.
    * @return True if the expression may be evaluated.
    */
  def validateOptionalExpression(wdlExpressionMaybe: Option[WdlValue]): Boolean = {
    wdlExpressionMaybe match {
      case None => staticDefaultOption.isDefined || validateNone.isValid
      case Some(wdlExpression: WdlExpression) =>
        /*
        TODO: BUG:

        Using `wdl4s.NoLookup` with the following options causes the following exception:

        command:
          sbt -J-Dbackend.default=SGE 'run run test_wdl/test.wdl - test_wdl/test.default.options'

        options:
          { "default_runtime_attributes": { "sge_queue": "fromoptions" } }

        exception:
          java.lang.RuntimeException: Expression evaluation failed due to java.lang.UnsupportedOperationException:
          No identifiers should be looked up: fromoptions: WdlExpression(<string:1:1 identifier "ZnJvbW9wdGlvbnM=">)

        Not sure why yet the string "fromoptions" is being converted to an expression and not a WdlString.

        For now, if something tries to "lookup" a value, convert it to a WdlString.
         */
        val wdlStringLookup: ScopedLookupFunction = (value: String) => WdlString(value)
        wdlExpression.evaluate(wdlStringLookup, PureStandardLibraryFunctions) match {
          case Success(wdlValue) => validateExpression.applyOrElse(wdlValue, (_: Any) => false)
          case Failure(throwable) =>
            throw new RuntimeException(s"Expression evaluation failed due to $throwable: $wdlExpression", throwable)
        }
      case Some(wdlValue) => validateExpression.applyOrElse(wdlValue, (_: Any) => false)
    }
  }

  /**
    * Returns an optional version of this validation.
    */
  final lazy val optional: OptionalRuntimeAttributesValidation[ValidatedType] =
  RuntimeAttributesValidation.optional(this)

  /**
    * Returns a version of this validation with the default value.
    *
    * @param wdlValue The default wdl value.
    * @return The new version of this validation.
    */
  final def withDefault(wdlValue: WdlValue) = RuntimeAttributesValidation.withDefault(this, wdlValue)

  /*
  Methods below provide aliases to expose protected methods to the package.
  Allows wrappers to wire their overrides to invoke the corresponding method on the inner object.
  The protected methods are only available to subclasses, or this trait. Now, no one outside this trait lineage can
  access the protected values, except the `validation` package that uses these back doors.
   */

  private[validation] lazy val validateValuePackagePrivate = validateValue

  private[validation] lazy val validateExpressionPackagePrivate = validateExpression

  private[validation] lazy val failureMessagePackagePrivate = failureMessage
}

/**
  * An optional version of a runtime attribute validation.
  *
  * @tparam ValidatedType The type of the validated value.
  */
trait OptionalRuntimeAttributesValidation[ValidatedType] extends RuntimeAttributesValidation[Option[ValidatedType]] {
  /**
    * Validates the wdl value.
    *
    * This method is the same as `validateValue`, but allows the implementor to not have to wrap the response in an
    * `Option`.
    *
    * @return The validated value or an error, wrapped in a cats validation.
    */
  protected def validateOption: PartialFunction[WdlValue, ErrorOr[ValidatedType]]

  override final protected lazy val validateValue = new PartialFunction[WdlValue, ErrorOr[Option[ValidatedType]]] {
    override def isDefinedAt(wdlValue: WdlValue) = validateOption.isDefinedAt(wdlValue)

    override def apply(wdlValue: WdlValue) = validateOption.apply(wdlValue).map(Option.apply)
  }

  override final protected lazy val validateNone = None.validNel[String]
}
