package cwl

import shapeless.Poly1
import wom.types.WomType

object RunOutputsToTypeMap extends Poly1 {

  def handleCommandLine(clt: CommandLineTool): Map[String, WomType] = {
    clt.outputs.toList.foldLeft(Map.empty[String, WomType]) {
      (acc, out) =>
        acc ++
          out.
            `type`.
            map(_.fold(MyriadOutputTypeToWomType)).
            map(out.id -> _).
            toList.
            toMap
    }
  }

  implicit def commandLineTool =
    at[CommandLineTool] {
      clt =>
          handleCommandLine(clt)
    }

  implicit def string = at[String] {
    fileName =>
      Map.empty[String, WomType]
  }

  implicit def expressionTool = at[ExpressionTool] {
    _ =>
        Map.empty[String, WomType]
  }

  implicit def workflow = at[Workflow] {
    wf =>  wf.steps.toList.flatMap(_.typedOutputs.toList).toMap
  }
}

