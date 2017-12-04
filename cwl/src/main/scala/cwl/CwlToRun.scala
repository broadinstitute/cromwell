package cwl

import shapeless._

object CwlToRun extends Poly1 {
  implicit def commandLineTool =
    at[CommandLineTool] {
      clt =>
       Coproduct[WorkflowStep.Run](clt)
    }

  implicit def workflow = at[Workflow] {
    wf =>
      Coproduct[WorkflowStep.Run](wf)
  }
}


