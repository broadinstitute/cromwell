package cwl

import shapeless.Poly1
import wom.types.WomType

object RunOutputsToTypeMap extends Poly1 {

  def handleOutputParameters[A <: OutputParameter](outputs: Array[A]): Map[String, WomType] = {
    outputs.toList.foldLeft(Map.empty[String, WomType]) {
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
        handleOutputParameters(clt.outputs)
    }

  implicit def string = at[String] {
    _ =>
      Map.empty[String, WomType]
  }

  implicit def expressionTool = at[ExpressionTool] {
    et =>
      handleOutputParameters(et.outputs)
  }

  implicit def workflow = at[Workflow] {
    wf =>
      handleOutputParameters(wf.outputs)
  }
}

