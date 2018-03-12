package cwl

import shapeless._
import WorkflowStep.Run

object RunToEmbeddedCwl extends Poly1 {
  implicit def commandLineTool =
    at[CommandLineTool] {
      clt =>
        (_: Map[String, Cwl]) =>
         Coproduct[Run](clt)
    }

  implicit def string = at[String] {
    fileName =>
      (cwlMap: Map[String, Cwl]) => {
        val cwl = cwlMap(fileName)
        cwl.fold(CwlToRun)
      }
  }


  implicit def expressionTool = at[ExpressionTool] {
    et =>
      (_: Map[String, Cwl]) =>
        Coproduct[Run](et)
  }

  implicit def workflow = at[Workflow] {
    wf =>
      (_: Map[String, Cwl]) =>
        Coproduct[Run](wf)
  }
}


