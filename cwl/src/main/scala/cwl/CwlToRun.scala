package cwl

import shapeless._

object CwlToRun extends Poly1 {
  implicit def commandLineTool = at[CommandLineTool] { Coproduct[WorkflowStep.Run](_) }
  implicit def workflow = at[Workflow] { Coproduct[WorkflowStep.Run](_) }
  implicit def expressionTool = at[ExpressionTool] { Coproduct[WorkflowStep.Run](_) }
}


