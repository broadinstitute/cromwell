package wdl4s.cwl

import org.scalatest.FlatSpec
/*
import shapeless.Coproduct
import shapeless.syntax.inject._
import shapeless.ops.coproduct.Inject
import mouse.all._
*/

class ThreeStepSpec extends FlatSpec {

  "Cwl Three step" should "convert to Wom" in {

    //val ps = new CommandLineTool()

    /*
    val patternInput =
      InputParameter(
        Some("pattern"),
        `type` = Some(Coproduct[MyriadInputType](CwlType.String)))
    val _inputs = Coproduct[WorkflowInput](Array(patternInput))

    val outputCgrep =
      WorkflowOutputParameter(
        id = Some("cgrep.count"),
        `type` = Some(Coproduct[MyriadOutputType](CwlType.Int)))

    val outputWc =
      WorkflowOutputParameter(
        id = Some("wc.count"),
        `type` = Some(Coproduct[MyriadOutputType](CwlType.Int)))

    //val _outputs = Coproduct[WorkflowInput](Array(outputCgrep, outputWc))

    val psWfStep  = WorkflowStep(
      id = Some("ps"),
      inputs = Coproduct[CommandLineTool.Inputs](Array.empty),
    )

    val m = new Workflow(
      None,
      `class` = "Workflow",
      inputs = _inputs,
      outputs = _outputs
      steps = Array.empty[WorkflowStep].inject[WorkflowSteps])
    */
  }


}
