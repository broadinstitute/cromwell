package cwl

import shapeless.Poly1

object RunToInputTypeMap extends Poly1 {

  type MyriadInputTypeMap = Map[String, Option[MyriadInputType]]

  import Case.Aux

  implicit def s: Aux[String, MyriadInputTypeMap] =
    at[String] {
      _ => throw new RuntimeException("this embedded CWL was supposed to have been embedded")
    }

  implicit def clt: Aux[CommandLineTool, MyriadInputTypeMap] =
    at[CommandLineTool]{_.inputs.map{
      input =>
        input.id -> input.`type`
    }.toMap}

  implicit def et: Aux[ExpressionTool, MyriadInputTypeMap] =
    at[ExpressionTool]{_.inputs.map{
      input =>
        input.id -> input.`type`
    }.toMap}

  implicit def wf: Aux[Workflow, MyriadInputTypeMap] =
    at[Workflow]{_.inputs.map{
      input =>
        input.id -> input.`type`
    }.toMap}

}
