package cromwell.backend

import _root_.wdl.draft2.model._
import cromwell.core.WorkflowOptions
import cromwell.util.JsonFormatting.WomValueJsonFormatter
import common.validation.ErrorOr.ErrorOr
import wom.callable.Callable.InputDefinition
import wom.expression.IoFunctionSet
import wom.values.WomObject
import wom.values.WomValue
import wom.{RuntimeAttributes, WomExpressionException}

import scala.util.{Success, Try}

/**
  * @param name Attribute name (LHS of name: "value" in the runtime section).
  * @param factoryDefault An optional default value for this attribute.
  * @param usedInCallCaching Whether or not this attribute is used to determine a cache hit or miss. ALL attributes are hashed and stored.
  */
case class RuntimeAttributeDefinition(name: String, factoryDefault: Option[WomValue], usedInCallCaching: Boolean)

object RuntimeAttributeDefinition {

  def evaluateRuntimeAttributes(unevaluated: RuntimeAttributes,
                                wdlFunctions: IoFunctionSet,
                                evaluatedInputs: Map[InputDefinition, WomValue],
                                platform: Option[Platform] = None
  ): ErrorOr[Map[String, WomValue]] = {
    import common.validation.ErrorOr._
    val inputsMap = evaluatedInputs map { case (x, y) => x.name -> y }
    val evaluated = unevaluated.attributes.traverseValues(_.evaluateValue(inputsMap, wdlFunctions))
    evaluated.map(e => applyPlatform(e, platform))
  }

  def applyPlatform(attributes: Map[String, WomValue], maybePlatform: Option[Platform]): Map[String, WomValue] = {

    def extractPlatformAttributes(platform: Platform): Map[String, WomValue] = {
      attributes.get(platform.runtimeKey) match {
        case Some(obj: WomObject) =>
          // WDL spec: "Use objects to avoid collisions"
          // https://github.com/openwdl/wdl/blob/wdl-1.1/SPEC.md#conventions-and-best-practices
          obj.values
        case _ =>
          // A malformed non-object override such as "gcp": "banana" is ignored
          Map.empty
      }
    }

    val platformAttributes = maybePlatform match {
      case Some(platform) =>
        extractPlatformAttributes(platform)
      case None =>
        Map.empty
    }

    // We've scooped our desired platform, now delete "azure", "gcp", etc.
    val originalAttributesWithoutPlatforms: Map[String, WomValue] =
      attributes -- Platform.all.map(_.runtimeKey)

    // With `++` keys from the RHS overwrite duplicates in LHS, which is what we want
    // RHS `Map.empty` is a no-op
    originalAttributesWithoutPlatforms ++ platformAttributes
  }

  def buildMapBasedLookup(evaluatedDeclarations: Map[InputDefinition, Try[WomValue]])(identifier: String): WomValue = {
    val successfulEvaluations = evaluatedDeclarations collect {
      case (k, v) if v.isSuccess => k.name -> v.get
    }
    successfulEvaluations.getOrElse(
      identifier,
      throw new WomExpressionException(s"Could not resolve variable $identifier as a task input")
    )
  }

  def addDefaultsToAttributes(runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
                              workflowOptions: WorkflowOptions
  )(specifiedAttributes: Map[LocallyQualifiedName, WomValue]): Map[LocallyQualifiedName, WomValue] = {
    import WomValueJsonFormatter._

    def isUnspecifiedAttribute(name: String) = !specifiedAttributes.contains(name)

    val missing = runtimeAttributeDefinitions filter { x => isUnspecifiedAttribute(x.name) }
    val defaults = missing map { x => (x, workflowOptions.getDefaultRuntimeOption(x.name)) } collect {
      case (runtimeAttributeDefinition, Success(jsValue)) =>
        runtimeAttributeDefinition.name -> jsValue.convertTo[WomValue]
      case (RuntimeAttributeDefinition(name, Some(factoryDefault), _), _) => name -> factoryDefault
    }
    specifiedAttributes ++ defaults
  }
}
