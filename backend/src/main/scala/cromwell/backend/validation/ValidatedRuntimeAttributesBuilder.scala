package cromwell.backend.validation

import cromwell.backend.RuntimeAttributeDefinition
import cromwell.core._
import lenthall.exception.MessageAggregation
import org.slf4j.Logger
import wdl4s.types.WdlType
import wdl4s.values.WdlValue

import scalaz.{Failure, Success}

final case class ValidatedRuntimeAttributes(attributes: Map[String, Any])

/**
  * Uses a collection of `RuntimeAttributesValidation` to build a `ValidatedRuntimeAttributes`, or throw an exception
  * with all of the validation errors.
  */
trait ValidatedRuntimeAttributesBuilder {

  /**
    * Returns the validations that should be used and returned during evaluation.
    *
    * These validations will be checked on backend initialization against all the calls. If a validation value will not
    * be supported, but should be checked anyway, supply it in `valuesOnlyValidations`, not here.
    *
    * @return the validations that should be used and returned during evaluation.
    */
  def validations: Seq[RuntimeAttributesValidation[_]]

  /**
    * Returns a mapping of the validations: RuntimeAttributesValidation each converted to a RuntimeAttributeDefinition.
    */
  final lazy val definitions: Seq[RuntimeAttributeDefinition] = {
    validations map RuntimeAttributesValidation.toRuntimeAttributeDefinition
  }

  /**
    * Returns the additional validations that should be used during value parsing.
    *
    * For example, sometimes docker might not be supported, BUT we want to still validate the value if specified.
    *
    * In that case, return the validation here.
    *
    * @return the additional validations that should be used during value parsing.
    */
  protected def unsupportedExtraValidations: Seq[OptionalRuntimeAttributesValidation[_]] = Seq.empty

  def unsupportedKeys(keys: Seq[String]): Seq[String] = keys.diff(validationKeys)

  private lazy val coercionMap: Map[String, Traversable[WdlType]] =
    validations.map(attr => attr.key -> attr.coercion).toMap

  private lazy val validationKeys = validations.map(_.key)

  private lazy val staticDefaults: Map[String, WdlValue] = (validations collect {
    case attr if attr.staticDefaultOption.isDefined => attr.key -> attr.staticDefaultOption.get
  }).toMap

  def build(attrs: Map[String, WdlValue], options: WorkflowOptions, logger: Logger): ValidatedRuntimeAttributes = {
    import RuntimeAttributesDefault._

    // Fail now if some workflow options are specified but can't be parsed correctly
    val defaultFromOptions = workflowOptionsDefault(options, coercionMap).get
    val withDefaultValues: Map[String, WdlValue] = withDefaults(attrs, List(defaultFromOptions, staticDefaults))

    RuntimeAttributesValidation.warnUnrecognized(withDefaultValues.keySet, validationKeys.toSet, logger)

    val runtimeAttributesErrorOr: ErrorOr[ValidatedRuntimeAttributes] = validate(withDefaultValues)
    runtimeAttributesErrorOr match {
      case Success(runtimeAttributes) => runtimeAttributes
      case Failure(nel) => throw new RuntimeException with MessageAggregation {
        override def exceptionContext: String = "Runtime attribute validation failed"

        override def errorMessages: Traversable[String] = nel.list
      }
    }
  }

  private def validate(values: Map[String, WdlValue]): ErrorOr[ValidatedRuntimeAttributes] = {
    val validationsForValues: Seq[RuntimeAttributesValidation[_]] = validations ++ unsupportedExtraValidations
    val errorsOrValuesMap: Seq[(String, ErrorOr[Any])] =
      validationsForValues.map(validation => validation.key -> validation.validate(values))

    import scalaz.Scalaz._

    val emptyResult: ErrorOr[List[(String, Any)]] = List.empty[(String, Any)].success
    val validationResult = errorsOrValuesMap.foldLeft(emptyResult) { (agg, errorOrValue) =>
      agg +++ {
        errorOrValue match {
          case (key, Success(value)) => List(key -> value).success
          case (key, Failure(nel)) => nel.failure
        }
      }
    }

    validationResult.map(result => ValidatedRuntimeAttributes(result.toMap))
  }
}
