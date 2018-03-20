package cromwell.backend.validation

import cats.data.Validated.{Invalid, Valid}
import cats.data.{NonEmptyList, Validated}
import cats.syntax.validated._
import com.typesafe.config.Config
import cromwell.backend.RuntimeAttributeDefinition
import common.validation.ErrorOr._
import org.slf4j.Logger
import wdl.draft2.model.expression.PureStandardLibraryFunctions
import wdl.draft2.model.{NoLookup, WdlExpression}
import wom.expression.{NoIoFunctionSet, WomExpression}
import wom.format.MemorySize
import wom.types._
import wom.values._

import scala.util.{Failure, Success}

object RuntimeAttributesValidation {

  def warnUnrecognized(actual: Set[String], expected: Set[String], logger: Logger): Unit = {
    val unrecognized = actual.diff(expected).mkString(", ")
    if (unrecognized.nonEmpty) logger.warn(s"Unrecognized runtime attribute keys: $unrecognized")
  }

  def validateDocker(docker: Option[WomValue], onMissingKey: => ErrorOr[Option[String]]): ErrorOr[Option[String]] = {
    validateWithValidation(docker, DockerValidation.instance.optional, onMissingKey)
  }

  def validateFailOnStderr(value: Option[WomValue], onMissingKey: => ErrorOr[Boolean]): ErrorOr[Boolean] = {
    validateWithValidation(value, FailOnStderrValidation.instance, onMissingKey)
  }

  def validateContinueOnReturnCode(value: Option[WomValue],
                                   onMissingKey: => ErrorOr[ContinueOnReturnCode]): ErrorOr[ContinueOnReturnCode] = {
    validateWithValidation(value, ContinueOnReturnCodeValidation.instance, onMissingKey)
  }

  def validateMemory(value: Option[WomValue], onMissingKey: => ErrorOr[MemorySize]): ErrorOr[MemorySize] = {
    validateWithValidation(value, MemoryValidation.instance(), onMissingKey)
  }

  def validateCpu(cpu: Option[WomValue], onMissingKey: => ErrorOr[Int]): ErrorOr[Int] = {
    validateWithValidation(cpu, CpuValidation.instance, onMissingKey)
  }

  private def validateWithValidation[T](valueOption: Option[WomValue],
                                        validation: RuntimeAttributesValidation[T],
                                        onMissingValue: => ErrorOr[T]): ErrorOr[T] = {
    valueOption match {
      case Some(value) =>
        validation.validateValue.applyOrElse(value, (_: Any) => validation.invalidValueFailure(value))
      case None => onMissingValue
    }
  }

  def validateInt(value: WomValue): ErrorOr[Int] = {
    WomIntegerType.coerceRawValue(value) match {
      case scala.util.Success(WomInteger(i)) => i.intValue.validNel
      case _ => s"Could not coerce ${value.valueString} into an integer".invalidNel
    }
  }

  def validateBoolean(value: WomValue): ErrorOr[Boolean] = {
    WomBooleanType.coerceRawValue(value) match {
      case scala.util.Success(WomBoolean(b)) => b.booleanValue.validNel
      case _ => s"Could not coerce ${value.valueString} into a boolean".invalidNel
    }
  }

  def parseMemoryString(k: String, s: WomString): ErrorOr[MemorySize] = {
    MemoryValidation.validateMemoryString(k, s)
  }

  def parseMemoryInteger(k: String, i: WomInteger): ErrorOr[MemorySize] = {
    MemoryValidation.validateMemoryInteger(k, i)
  }

  def withDefault[ValidatedType](validation: RuntimeAttributesValidation[ValidatedType],
                                 default: WomValue): RuntimeAttributesValidation[ValidatedType] = {
    new RuntimeAttributesValidation[ValidatedType] {
      override def key: String = validation.key

      override def coercion: Traversable[WomType] = validation.coercion

      override protected def validateValue: PartialFunction[WomValue, ErrorOr[ValidatedType]] =
        validation.validateValuePackagePrivate

      override protected def validateExpression: PartialFunction[WomValue, Boolean] =
        validation.validateExpressionPackagePrivate

      override protected def invalidValueMessage(value: WomValue): String =
        validation.invalidValueMessagePackagePrivate(value)

      override protected def missingValueMessage: String = validation.missingValueMessage

      override protected def usedInCallCaching: Boolean = validation.usedInCallCachingPackagePrivate

      override protected def staticDefaultOption = Option(default)
    }
  }

  def optional[ValidatedType](validation: RuntimeAttributesValidation[ValidatedType]):
  OptionalRuntimeAttributesValidation[ValidatedType] = {
    new OptionalRuntimeAttributesValidation[ValidatedType] {
      override def key: String = validation.key

      override def coercion: Traversable[WomType] = validation.coercion

      override protected def validateOption: PartialFunction[WomValue, ErrorOr[ValidatedType]] =
        validation.validateValuePackagePrivate

      override protected def validateExpression: PartialFunction[WomValue, Boolean] =
        validation.validateExpressionPackagePrivate

      override protected def invalidValueMessage(value: WomValue): String =
        validation.invalidValueMessagePackagePrivate(value)

      override protected def missingValueMessage: String = validation.missingValueMessage

      override protected def usedInCallCaching: Boolean = validation.usedInCallCachingPackagePrivate

      override protected def staticDefaultOption = validation.staticDefaultOption
    }
  }

  /**
    * Returns the value from the attributes, unpacking options, and converting them to string values suitable for
    * storage in metadata.
    *
    * @param validatedRuntimeAttributes The values to search.
    * @return The keys and extracted values.
    */
  def toMetadataStrings(validatedRuntimeAttributes: ValidatedRuntimeAttributes): Map[String, String] = {
    val attributeOptions: Map[String, Option[Any]] = validatedRuntimeAttributes.attributes.mapValues(unpackOption)

    val attributes: Map[String, String] = attributeOptions collect {
      case (name, Some(values: Traversable[_])) => (name, values.mkString(","))
      case (name, Some(value)) => (name, value.toString)
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
}

/**
  * A wrapper class to classify config-based default runtime attributes
  * that cannot be coerced into an acceptable WomType.
  */
case class BadDefaultAttribute(badDefaultValue: WomValue) extends WomValue {
  val womType = WomStringType
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
  def coercion: Traversable[WomType]

  /**
    * Validates the wdl value.
    *
    * @return The validated value or an error, wrapped in a cats validation.
    */
  protected def validateValue: PartialFunction[WomValue, ErrorOr[ValidatedType]]

  /**
    * Returns the value for when there is no wdl value. By default returns an error.
    *
    * @return the value for when there is no wdl value.
    */
  protected def validateNone: ErrorOr[ValidatedType] = missingValueFailure

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
  protected def validateExpression: PartialFunction[WomValue, Boolean] = {
    case womValue => coercion.exists(_ == womValue.womType)
  }

  /**
    * Returns the optional default value when no other is specified.
    *
    * @return the optional default value when no other is specified.
    */
  protected def staticDefaultOption: Option[WomValue] = None

  /**
    * Returns message to return when a value is invalid.
    *
    * By default returns the missingValueMessage.
    *
    * @return Message to return when a value is invalid.
    */
  protected def invalidValueMessage(value: WomValue): String = missingValueMessage

  /**
    * Utility method to wrap the invalidValueMessage in an ErrorOr.
    *
    * @return Wrapped invalidValueMessage.
    */
  protected final def invalidValueFailure(value: WomValue): ErrorOr[ValidatedType] =
    invalidValueMessage(value).invalidNel

  /**
    * Returns message to return when a value is missing.
    *
    * @return Message to return when a value is missing.
    */
  protected def missingValueMessage: String = s"Expecting $key runtime attribute to be a type in $coercion"

  /**
    * Utility method to wrap the missingValueMessage in an ErrorOr.
    *
    * @return Wrapped missingValueMessage.
    */
  protected final lazy val missingValueFailure: ErrorOr[ValidatedType] = missingValueMessage.invalidNel

  /**
    * Runs this validation on the value matching key.
    *
    * NOTE: The values passed to this method should already be evaluated instances of WomValue, and not WomExpression.
    *
    * @param values The full set of values.
    * @return The error or valid value for this key.
    */
  def validate(values: Map[String, WomValue]): ErrorOr[ValidatedType] = {
    values.get(key) match {
      case Some(value) => validateValue.applyOrElse(value, (_: Any) => invalidValueFailure(value))
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
    * NOTE: If there is an attempt lookup a value within a WdlExpression, or a WdlExpression fails to evaluate for any
    * reason, this method will simply return true.
    *
    * @param wdlExpressionMaybe The optional expression.
    * @return True if the expression may be evaluated.
    */
  def validateOptionalWomValue(wdlExpressionMaybe: Option[WomValue]): Boolean = {
    wdlExpressionMaybe match {
      case None => staticDefaultOption.isDefined || validateNone.isValid
      case Some(wdlExpression: WdlExpression) =>
        wdlExpression.evaluate(NoLookup, PureStandardLibraryFunctions) match {
          case Success(womValue) => validateExpression.applyOrElse(womValue, (_: Any) => false)
          case Failure(_) => true // If we can't evaluate it, we'll let it pass for now...
        }
      case Some(womValue) => validateExpression.applyOrElse(womValue, (_: Any) => false)
    }
  }

  def validateOptionalWomExpression(womExpressionMaybe: Option[WomExpression]): Boolean = {
    womExpressionMaybe match {
      case None => staticDefaultOption.isDefined || validateNone.isValid
      case Some(womExpression) =>
        womExpression.evaluateValue(Map.empty, NoIoFunctionSet) match {
          case Valid(womValue) => validateExpression.applyOrElse(womValue, (_: Any) => false)
          case Invalid(_) => true // If we can't evaluate it, we'll let it pass for now...
        }
    }
  }
  /**
    * Used to convert this instance to a `RuntimeAttributeDefinition`.
    *
    * @see [[RuntimeAttributeDefinition.usedInCallCaching]]
    * @return Value for [[RuntimeAttributeDefinition.usedInCallCaching]].
    */
  protected def usedInCallCaching: Boolean = false

  /**
    * Returns this as an instance of a runtime attribute definition.
    */
  final lazy val runtimeAttributeDefinition = RuntimeAttributeDefinition(key, staticDefaultOption, usedInCallCaching)

  /**
    * Returns an optional version of this validation.
    */
  final lazy val optional: OptionalRuntimeAttributesValidation[ValidatedType] =
  RuntimeAttributesValidation.optional(this)

  /**
    * Returns a version of this validation with the default value.
    *
    * @param womValue The default wdl value.
    * @return The new version of this validation.
    */
  final def withDefault(womValue: WomValue): RuntimeAttributesValidation[ValidatedType] =
    RuntimeAttributesValidation.withDefault(this, womValue)

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
  final def configDefaultWomValue(optionalRuntimeConfig: Option[Config]): Option[WomValue] = {
    optionalRuntimeConfig flatMap { config =>
      val value = config.getValue(key).unwrapped()
      coercion.collectFirst({
        case womType if womType.coerceRawValue(value).isSuccess => {
          womType.coerceRawValue(value).get
        }
      }) orElse Option(BadDefaultAttribute(WomString(value.toString)))
    }
  }

  final def configDefaultValue(optionalRuntimeConfig: Option[Config]): Option[String] = {
    optionalRuntimeConfig match {
      case Some(config) if config.hasPath(key) => Option(config.getValue(key).unwrapped().toString)
      case _ => None
    }
  }

  /*
  Methods below provide aliases to expose protected methods to the package.
  Allows wrappers to wire their overrides to invoke the corresponding method on the inner object.
  The protected methods are only available to subclasses, or this trait. Now, no one outside this trait lineage can
  access the protected values, except the `validation` package that uses these back doors.
   */

  private[validation] final lazy val validateValuePackagePrivate = validateValue

  private[validation] final lazy val validateExpressionPackagePrivate = validateExpression

  private[validation] final def invalidValueMessagePackagePrivate(value: WomValue) = invalidValueMessage(value)

  private[validation] final lazy val missingValueMessagePackagePrivate = missingValueMessage

  private[validation] final lazy val usedInCallCachingPackagePrivate = usedInCallCaching
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
  protected def validateOption: PartialFunction[WomValue, ErrorOr[ValidatedType]]

  override final protected lazy val validateValue = new PartialFunction[WomValue, ErrorOr[Option[ValidatedType]]] {
    override def isDefinedAt(womValue: WomValue): Boolean = validateOption.isDefinedAt(womValue)

    override def apply(womValue: WomValue): Validated[NonEmptyList[String], Option[ValidatedType]] = {
      validateOption.apply(womValue).map(Option.apply)
    }
  }

  override final protected lazy val validateNone: ErrorOr[None.type] = None.validNel[String]
}
