package cwl

import shapeless.Poly1
import cwl.CwlType.CwlType
import wom.types.WomType

object RunOutputsToTypeMap extends Poly1 {

  def handleCommandLine(clt: CommandLineTool): Map[String, WomType] = {
    clt.outputs.toList.foldLeft(Map.empty[String, WomType]) {
      (acc, out) =>
        acc ++
          out.
            `type`.
            flatMap(_.select[CwlType]).
            map(cwlTypeToWomType).
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

