package cwl

import cwl.command.ParentName
import shapeless.Poly1

object RunToInputTypeMap extends Poly1 {

  type MyriadInputTypeMap = Map[String, Option[MyriadInputType]]

  type OutputType = ParentName => MyriadInputTypeMap

  import Case.Aux

  implicit def s: Aux[String, OutputType] =
    at[String] {
      run =>  _ => throw new RuntimeException(s"Run field $run was not inlined as expected")
    }

  implicit def clt: Aux[CommandLineTool, OutputType] =
    at[CommandLineTool]{clt => parentName => clt.inputs.map{
      input =>
          FullyQualifiedName(input.id)(parentName).id -> input.`type`
    }.toMap}

  implicit def et: Aux[ExpressionTool, OutputType] =
    at[ExpressionTool]{et => parentName => et.inputs.map{
      input =>
        FullyQualifiedName(input.id)(parentName).id -> input.`type`
    }.toMap}

  implicit def wf: Aux[Workflow, OutputType] =
    at[Workflow]{ wf => parentName => wf.inputs.map{
      input =>
        FullyQualifiedName(input.id)(parentName).id -> input.`type`
    }.toMap}

}
