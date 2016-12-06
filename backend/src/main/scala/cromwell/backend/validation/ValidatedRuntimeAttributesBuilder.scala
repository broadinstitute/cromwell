package cromwell.backend.validation

import cats.data.Validated._
import cats.instances.list._
import cromwell.backend.RuntimeAttributeDefinition
import lenthall.exception.MessageAggregation
import lenthall.validation.ErrorOr._
import org.slf4j.Logger
import wdl4s.values.WdlValue

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

  private lazy val validationKeys = validations.map(_.key)

  def build(attrs: Map[String, WdlValue], logger: Logger): ValidatedRuntimeAttributes = {
    RuntimeAttributesValidation.warnUnrecognized(attrs.keySet, validationKeys.toSet, logger)

    val runtimeAttributesErrorOr: ErrorOr[ValidatedRuntimeAttributes] = validate(attrs)
    runtimeAttributesErrorOr match {
      case Valid(runtimeAttributes) => runtimeAttributes
      case Invalid(nel) => throw new RuntimeException with MessageAggregation {
        override def exceptionContext: String = "Runtime attribute validation failed"

        override def errorMessages: Traversable[String] = nel.toList
      }
    }
  }

  private def validate(values: Map[String, WdlValue]): ErrorOr[ValidatedRuntimeAttributes] = {
    val validationsForValues = validations ++ unsupportedExtraValidations
    val listOfKeysToErrorOrAnys: List[(String, ErrorOr[Any])] =
      validationsForValues.map(validation => validation.key -> validation.validate(values)).toList

    val listOfErrorOrKeysToAnys: List[ErrorOr[(String, Any)]] = listOfKeysToErrorOrAnys map {
      case (key, errorOrAny) => errorOrAny map { any => (key, any) }
    }

    import cats.syntax.traverse._
    val errorOrListOfKeysToAnys: ErrorOr[List[(String, Any)]] = listOfErrorOrKeysToAnys.sequence[ErrorOr, (String, Any)]
    errorOrListOfKeysToAnys map { listOfKeysToAnys => ValidatedRuntimeAttributes(listOfKeysToAnys.toMap) }
  }
}
