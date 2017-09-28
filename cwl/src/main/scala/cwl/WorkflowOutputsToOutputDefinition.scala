package cwl

import shapeless.Poly1
import wom.callable.Callable.OutputDefinition
import wom.expression.PlaceholderWomExpression

object WorkflowOutputsToOutputDefinition extends Poly1 {

  def fullIdToOutputDefintition(fullyQualifiedName: String, typeMap: WdlTypeMap) = {

    //we want to only look at the id, not the filename
    val lookupId = WorkflowStepInputOrOutputId(fullyQualifiedName).ioId

    OutputDefinition(fullyQualifiedName, typeMap(lookupId), PlaceholderWomExpression(Set.empty, typeMap(lookupId)))
  }

  implicit def a = at[Array[WorkflowStepOutput]] { outputs =>
    (typeMap: WdlTypeMap) =>
      outputs.map(output => fullIdToOutputDefintition(output.id, typeMap)).toSet
  }

  implicit def b = at[Array[String]] { outputs =>
    (typeMap: WdlTypeMap) =>
      outputs.map(fullIdToOutputDefintition(_, typeMap)).toSet
  }

}

