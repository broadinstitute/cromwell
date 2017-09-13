package cromwell.backend

import cromwell.core.{NoIoFunctionSet, WorkflowOptions}
import cromwell.util.JsonFormatting.WdlValueJsonFormatter
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wdl.values.WdlValue
import wdl4s.wdl.{WdlExpressionException, _}
import wdl4s.wom.RuntimeAttributes
import wdl4s.wom.callable.Callable.InputDefinition
import wdl4s.wom.expression.IoFunctionSet

import scala.util.{Success, Try}

/**
  * @param name Attribute name (LHS of name: "value" in the runtime section).
  * @param factoryDefault An optional default value for this attribute.
  * @param usedInCallCaching Whether or not this attribute is used to determine a cache hit or miss. ALL attributes are hashed and stored.
  */
case class RuntimeAttributeDefinition(name: String, factoryDefault: Option[WdlValue], usedInCallCaching: Boolean)

object RuntimeAttributeDefinition {

  def evaluateRuntimeAttributes(unevaluated: RuntimeAttributes,
                                wdlFunctions: IoFunctionSet,
                                evaluatedInputs: Map[InputDefinition, WdlValue]): ErrorOr[Map[String, WdlValue]] = {
    import lenthall.validation.ErrorOr._
    val inputsMap = evaluatedInputs map { case (x, y) => x.name -> y }
    unevaluated.attributes
      .map({ case (key, expression) => key -> expression.evaluateValue(inputsMap, NoIoFunctionSet)})
      .sequence
  }

  def buildMapBasedLookup(evaluatedDeclarations: Map[InputDefinition, Try[WdlValue]])(identifier: String): WdlValue = {
    val successfulEvaluations = evaluatedDeclarations collect {
      case (k, v) if v.isSuccess => k.name -> v.get
    }
    successfulEvaluations.getOrElse(identifier, throw new WdlExpressionException(s"Could not resolve variable $identifier as a task input"))
  }

  def addDefaultsToAttributes(runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition], workflowOptions: WorkflowOptions)
                             (specifiedAttributes: Map[LocallyQualifiedName, WdlValue]): Map[LocallyQualifiedName, WdlValue] = {
    import WdlValueJsonFormatter._

    def isUnspecifiedAttribute(name: String) = !specifiedAttributes.contains(name)

    val missing = runtimeAttributeDefinitions filter { x => isUnspecifiedAttribute(x.name) }
    val defaults = missing map { x => (x, workflowOptions.getDefaultRuntimeOption(x.name)) } collect {
      case (runtimeAttributeDefinition, Success(jsValue)) => runtimeAttributeDefinition.name -> jsValue.convertTo[WdlValue]
      case (RuntimeAttributeDefinition(name, Some(factoryDefault), _), _) => name -> factoryDefault
    }
    specifiedAttributes ++ defaults
  }
}
