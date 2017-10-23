package cromwell.backend

import _root_.wdl._
import cromwell.core.{NoIoFunctionSet, WorkflowOptions}
import cromwell.util.JsonFormatting.WdlValueJsonFormatter
import lenthall.validation.ErrorOr.ErrorOr
import wom.callable.Callable.InputDefinition
import wom.expression.IoFunctionSet
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
                                evaluatedInputs: Map[InputDefinition, WomValue]): ErrorOr[Map[String, WomValue]] = {
    import lenthall.validation.ErrorOr._
    val inputsMap = evaluatedInputs map { case (x, y) => x.name -> y }
    unevaluated.attributes.traverseValues(_.evaluateValue(inputsMap, NoIoFunctionSet))
  }

  def buildMapBasedLookup(evaluatedDeclarations: Map[InputDefinition, Try[WomValue]])(identifier: String): WomValue = {
    val successfulEvaluations = evaluatedDeclarations collect {
      case (k, v) if v.isSuccess => k.name -> v.get
    }
    successfulEvaluations.getOrElse(identifier, throw new WomExpressionException(s"Could not resolve variable $identifier as a task input"))
  }

  def addDefaultsToAttributes(runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition], workflowOptions: WorkflowOptions)
                             (specifiedAttributes: Map[LocallyQualifiedName, WomValue]): Map[LocallyQualifiedName, WomValue] = {
    import WdlValueJsonFormatter._

    def isUnspecifiedAttribute(name: String) = !specifiedAttributes.contains(name)

    val missing = runtimeAttributeDefinitions filter { x => isUnspecifiedAttribute(x.name) }
    val defaults = missing map { x => (x, workflowOptions.getDefaultRuntimeOption(x.name)) } collect {
      case (runtimeAttributeDefinition, Success(jsValue)) => runtimeAttributeDefinition.name -> jsValue.convertTo[WomValue]
      case (RuntimeAttributeDefinition(name, Some(factoryDefault), _), _) => name -> factoryDefault
    }
    specifiedAttributes ++ defaults
  }
}
