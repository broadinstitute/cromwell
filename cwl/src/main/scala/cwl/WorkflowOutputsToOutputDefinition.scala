package cwl

import shapeless.Poly1
import wom.callable.Callable.OutputDefinition
import wom.expression.WomExpression

object WorkflowOutputsToOutputDefinition extends Poly1 {

  def fullIdToOutputDefintition(fullyQualifiedName: String, typeMap: WdlTypeMap, expression: Map[String, WomExpression]) = {

    //we want to only look at the id, not the filename
    val lookupId = FullyQualifiedName(fullyQualifiedName).id

    OutputDefinition(fullyQualifiedName, typeMap(lookupId), expression(lookupId))
  }

  implicit def a = at[Array[WorkflowStepOutput]] { outputs =>
    (typeMap: WdlTypeMap, expressionMap: Map[String, WomExpression]) =>
      outputs.map(output => fullIdToOutputDefintition(output.id, typeMap, expressionMap)).toSet
  }

  implicit def b = at[Array[String]] { outputs =>
    (typeMap: WdlTypeMap, expressionMap: Map[String, WomExpression]) =>
      outputs.map(fullIdToOutputDefintition(_, typeMap, expressionMap)).toSet
  }

}

