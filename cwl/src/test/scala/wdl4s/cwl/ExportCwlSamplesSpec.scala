package wdl4s.cwl

import io.circe.yaml.Printer
import org.scalatest.{FlatSpec, Matchers}
import shapeless.{Coproduct, Witness}
import wdl4s.cwl.CommandLineTool.{BaseCommand, Outputs}
import wdl4s.cwl.CwlVersion.CwlVersion
import wdl4s.cwl.WorkflowStep.{Inputs, Run}

class ExportCwlSamplesSpec extends FlatSpec with Matchers {

  def assertCorrectJson(cwl: Cwl, expectedYaml: String): Unit = CwlCodecs.cwlToYaml(cwl) shouldBe expectedYaml

  it should "encode sample CWL command line tool" in {
    val tool =
      CommandLineTool(
        inputs = Coproduct[CommandLineTool.Inputs](Map("message" -> CommandInputParameter(
          inputBinding = Option(CommandLineBinding(
            position = Option(1)
          )),
          default = None,
          `type` = None
        ))),
        outputs = Coproduct[CommandLineTool.Outputs](Array.empty[CommandOutputParameter]),
        baseCommand = Option(Coproduct[BaseCommand]("echo"))
      )
    val expectedToolJsonString =
      """inputs:
        |  message:
        |    inputBinding:
        |      position: 1
        |outputs: []
        |class: CommandLineTool
        |cwlVersion: v1.0
        |baseCommand: echo
        |""".stripMargin
    assertCorrectJson(tool, expectedToolJsonString)
  }

  it should "encode sample CWL workflow" in {
    val workflow = Workflow(
      inputs = Coproduct[WorkflowInput](
        Map(
          "inp" -> Coproduct[MyriadInputType](CwlType.File),
          "ex" -> Coproduct[MyriadInputType](CwlType.String)
        )
      ),
      outputs = Coproduct[WorkflowOutput](
        Map(
          "classout" -> WorkflowOutputParameter(
            `type` = Option(Coproduct[MyriadOutputType](CwlType.File)),
            outputSource = Option(Coproduct[WorkflowOutputParameter#OutputSource]("compile/classfile"))
          )
        )
      ),
      steps = Coproduct[WorkflowSteps](
        Map(
          "untar" -> WorkflowStep(
            run = Coproduct[WorkflowStep.Run]("tar-param.cwl"),
            in = Coproduct[WorkflowStep.Inputs](
              Map(
                "tarfile" -> Coproduct[WorkflowStepInputSource]("inp"),
                "extractfile" -> Coproduct[WorkflowStepInputSource]("ex")
              )
            ),
            out = Coproduct[WorkflowStep.Outputs](Array("example_out"))
          ),
          "compile" -> WorkflowStep(
            run = Coproduct[WorkflowStep.Run]("arguments.cwl"),
            in = Coproduct[WorkflowStep.Inputs](
              Map(
                "src" -> Coproduct[WorkflowStepInputSource]("untar/example_out")
              )
            ),
            out = Coproduct[WorkflowStep.Outputs](Array("classfile"))
          )
        )
      ))
    val expectedWorkflowJsonString =
      """cwlVersion: v1.0
        |class: Workflow
        |inputs:
        |  inp: File
        |  ex: string
        |outputs:
        |  classout:
        |    outputSource: compile/classfile
        |    type: File
        |steps:
        |  untar:
        |    in:
        |      tarfile: inp
        |      extractfile: ex
        |    out:
        |    - example_out
        |    run: tar-param.cwl
        |  compile:
        |    in:
        |      src: untar/example_out
        |    out:
        |    - classfile
        |    run: arguments.cwl
        |""".stripMargin
    assertCorrectJson(workflow, expectedWorkflowJsonString)
  }

  it should "encode sample CWL env" in {
    val tool = CommandLineTool(
      baseCommand = Option(Coproduct[BaseCommand]("env")),
      requirements = Option(Array(Coproduct[Requirement](EnvVarRequirement(
        envDef = Coproduct[EnvVarRequirement.EnvDef](Map("HELLO" -> "$(inputs.message)"))
      )))) : Option[Array[Requirement]],
      inputs = Coproduct[CommandLineTool.Inputs](Map(
        "message" -> Coproduct[MyriadCommandInputType](CwlType.String)
      )),
      outputs = Coproduct[CommandLineTool.Outputs](Array[CommandOutputParameter]()) : CommandLineTool.Outputs
    )
    val expectedToolJsonString =
      """inputs:
        |  message: string
        |outputs: []
        |class: CommandLineTool
        |requirements:
        |- class: EnvVarRequirement
        |  envDef:
        |    HELLO: $(inputs.message)
        |cwlVersion: v1.0
        |baseCommand: env
        |""".stripMargin
    assertCorrectJson(tool, expectedToolJsonString)
  }

}
