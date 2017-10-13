package cwl

import shapeless.Poly1
import cwl.CwlType.CwlType
import wom.types.WdlType

object RunOutputsToTypeMap extends Poly1 {
  def mungeId(fullyQualifiedId: String): String = {
    val step = fullyQualifiedId.substring(fullyQualifiedId.lastIndexOf("#") + 1)
    // Doesn't matter if the string actually contains '/' or not, this takes the whole string if it's absent
    // which works out to be the right thing to do.
    step.substring(step.lastIndexOf("/") + 1)
  }

  def handleCommandLine(clt: CommandLineTool): Map[String, WdlType] = {
    clt.outputs.toList.foldLeft(Map.empty[String, WdlType]) {
      (acc, out) =>
        acc ++
          out.
            `type`.
            flatMap(_.select[CwlType]).
            map(cwlTypeToWdlType).
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
      Map.empty[String, WdlType]
  }

  implicit def expressionTool = at[ExpressionTool] {
    _ =>
        Map.empty[String, WdlType]
  }

  implicit def workflow = at[Workflow] {
    wf =>  wf.steps.toList.flatMap(_.typedOutputs.toList).toMap
  }
}

