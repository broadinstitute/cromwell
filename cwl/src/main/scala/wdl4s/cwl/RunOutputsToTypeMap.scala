package wdl4s.cwl

import shapeless.Poly1
import wdl4s.cwl.CwlType.CwlType
import wdl4s.wdl.types.WdlType

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
        (_: Map[String, CwlFile]) =>
          handleCommandLine(clt)
    }

  implicit def string = at[String] {
    fileName =>
      (cwlMap: Map[String, CwlFile]) =>
        cwlMap(fileName) match {
          case clt: CommandLineTool => handleCommandLine(clt)
          case wf: Workflow => handleWorkflow(wf)
        }
  }

  implicit def expressionTool = at[ExpressionTool] {
    _ =>
      (_: Map[String, CwlFile]) =>
        Map.empty[String, WdlType]
  }

  //TODO: run _.stepOutputs recursively
  def handleWorkflow(workflow: Workflow) = {
    Map.empty[String, WdlType]
  }

  implicit def workflow = at[Workflow] {
    wf =>
      (_: Map[String, CwlFile]) =>
        handleWorkflow(wf)
  }
}

