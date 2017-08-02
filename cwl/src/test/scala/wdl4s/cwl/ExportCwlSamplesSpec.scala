package wdl4s.cwl

import org.scalatest.{FlatSpec, Matchers}
import shapeless.Coproduct
import wdl4s.cwl.CommandLineTool.{BaseCommand, StringOrExpression}

class ExportCwlSamplesSpec extends FlatSpec with Matchers {

  def assertCorrectJson(cwl: Cwl, expectedYaml: String) = CwlCodecs.cwlToYaml(cwl) shouldBe expectedYaml

  it should "encode sample CWL command line tool" in {
    val tool =
      CommandLineTool(
        inputs =  Array(CommandInputParameter(
          id = "message",
          inputBinding = Option(CommandLineBinding(
            position = Option(1)
          )),
        )),
        baseCommand = Option(Coproduct[BaseCommand]("echo"))
      )
    val expectedToolJsonString =
"""inputs:
- id: message
  inputBinding:
    position: 1
outputs: []
class: CommandLineTool
cwlVersion: v1.0
baseCommand: echo
"""
    assertCorrectJson(tool, expectedToolJsonString)
  }

  it should "encode sample CWL workflow" in {
    val workflow = Workflow(
      inputs = Array(
          InputParameter(id = "inp", `type` = Option(Coproduct[MyriadInputType](CwlType.File))),
          InputParameter(id = "ex", `type` = Option(Coproduct[MyriadInputType](CwlType.String))),
        ),
      outputs = (
        Array(
           WorkflowOutputParameter(
            id = "classout",
            `type` = Option(Coproduct[MyriadOutputType](CwlType.File)),
            outputSource = Option(Coproduct[WorkflowOutputParameter#OutputSource]("compile/classfile"))
          )
        )
      ),
      steps =
        Array(
           WorkflowStep(
            id = "untar",
            run = Coproduct[WorkflowStep.Run]("tar-param.cwl"),
            in =
              Array(
                WorkflowStepInput(id = "tarfile", source = Option(Coproduct[WorkflowStepInputSource]("inp"))),
                WorkflowStepInput(id = "extractfile", source =  Option(Coproduct[WorkflowStepInputSource]("ex")))
              ),
            out = Coproduct[WorkflowStep.Outputs](Array("example_out"))
          ),
          WorkflowStep(
            id = "compile",
            run = Coproduct[WorkflowStep.Run]("arguments.cwl"),
            in = Array(
                 WorkflowStepInput(id = "src", source = Option(Coproduct[WorkflowStepInputSource]("untar/example_out")))
              ),
            out = Coproduct[WorkflowStep.Outputs](Array("classfile"))
          )
        )
      )
    val expectedWorkflowJsonString =
      """cwlVersion: v1.0
class: Workflow
inputs:
- id: inp
  type: File
- id: ex
  type: string
outputs:
- id: classout
  outputSource: compile/classfile
  type: File
steps:
- id: untar
  in:
  - id: tarfile
    source: inp
  - id: extractfile
    source: ex
  out:
  - example_out
  run: tar-param.cwl
- id: compile
  in:
  - id: src
    source: untar/example_out
  out:
  - classfile
  run: arguments.cwl
""".stripMargin
    assertCorrectJson(workflow, expectedWorkflowJsonString)
  }

  it should "encode sample CWL env" in {
    val tool = CommandLineTool(
      baseCommand = Option(Coproduct[BaseCommand]("env")),
      requirements = Option(Array(Coproduct[Requirement](EnvVarRequirement(
        envDef = Array(EnvironmentDef("HELLO", Coproduct[StringOrExpression]("$(inputs.message)")))
      )))),
      inputs = Array(
        CommandInputParameter(id =  "message", `type` = Option(Coproduct[MyriadInputType](CwlType.String)))
      )
    )
    val expectedToolJsonString =
"""inputs:
- id: message
  type: string
outputs: []
class: CommandLineTool
requirements:
- class: EnvVarRequirement
  envDef:
  - envName: HELLO
    envValue: $(inputs.message)
cwlVersion: v1.0
baseCommand: env
""".stripMargin
    assertCorrectJson(tool, expectedToolJsonString)
  }

}
