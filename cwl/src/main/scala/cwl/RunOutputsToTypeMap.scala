package cwl

import shapeless.Poly1
import wom.types.WomType

object RunOutputsToTypeMap extends Poly1 {

  import Case.Aux

  type SchemaDefToTypeMap = SchemaDefRequirement => WomTypeMap

  def handleOutputParameters[A <: OutputParameter](outputs: Array[A], schemaDefRequirement: SchemaDefRequirement): Map[String, WomType] = {
    outputs.toList.foldLeft(Map.empty[String, WomType]) {
      (acc, out) =>
        acc ++
          out.
            `type`.
            map(_.fold(MyriadOutputTypeToWomType).apply(schemaDefRequirement)).
            map(out.id -> _).
            toList.
            toMap
    }
  }

  implicit def commandLineTool : Aux[CommandLineTool, SchemaDefToTypeMap] =
    at[CommandLineTool] {
      clt => handleOutputParameters(clt.outputs, _)
    }

  implicit def string: Aux[String, SchemaDefToTypeMap] = at[String] {
    _ => _ =>
      Map.empty[String, WomType] //should never happen because we dereference/embed CWL subworkflows into one object from the original string references
  }

  implicit def expressionTool: Aux[ExpressionTool, SchemaDefToTypeMap] = at[ExpressionTool] {
    et =>
      handleOutputParameters(et.outputs, _)
  }

  implicit def workflow: Aux[Workflow, SchemaDefToTypeMap] = at[Workflow] {
    wf =>
      handleOutputParameters(wf.outputs, _)
  }
}

