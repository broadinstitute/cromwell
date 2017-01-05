package cromwell.backend

import cromwell.core.WorkflowOptions
import cromwell.util.JsonFormatting.WdlValueJsonFormatter
import lenthall.util.TryUtil
import wdl4s.{WdlExpressionException, _}
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values.WdlValue

import scala.util.{Success, Try}

/**
  * @param name Attribute name (LHS of name: "value" in the runtime section).
  * @param factoryDefault An optional default value for this attribute.
  * @param usedInCallCaching Whether or not this attribute is used to determine a cache hit or miss. ALL attributes are hashed and stored.
  */
case class RuntimeAttributeDefinition(name: String, factoryDefault: Option[WdlValue], usedInCallCaching: Boolean)

object RuntimeAttributeDefinition {

  def evaluateRuntimeAttributes(unevaluated: RuntimeAttributes,
                                wdlFunctions: WdlStandardLibraryFunctions,
                                evaluatedInputs: Map[Declaration, WdlValue]): Try[Map[String, WdlValue]] = {
    val tryInputs = evaluatedInputs map { case (x, y) => x.unqualifiedName -> Success(y) }
    val mapBasedLookup = buildMapBasedLookup(tryInputs) _
    val mapOfTries = unevaluated.attrs mapValues {
      expr => expr.evaluate(mapBasedLookup, wdlFunctions)
    }
    TryUtil.sequenceMap(mapOfTries)
  }

  def buildMapBasedLookup(evaluatedDeclarations: Map[LocallyQualifiedName, Try[WdlValue]])(identifier: String): WdlValue = {
    val successfulEvaluations = evaluatedDeclarations collect {
      case (k, v) if v.isSuccess => k -> v.get
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
