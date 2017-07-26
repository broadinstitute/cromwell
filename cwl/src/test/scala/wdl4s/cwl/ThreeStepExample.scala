package wdl4s.cwl

import shapeless._
import syntax.singleton._
import wdl4s.cwl.CommandLineTool.{Argument, BaseCommand, Inputs, StringOrExpression}
import wdl4s.cwl.CommandOutputBinding.Glob
import wdl4s.cwl.WorkflowStep.{Outputs, Run}
import io.circe.syntax._
import io.circe.yaml._
import io.circe.yaml.syntax._
import CwlCodecs._

/*
 This example calls `ps` , then counts the number of processes that match a pattern input.
 The output of the ps call is also counted via `wc`.
 */
object ThreeStepExample extends App {
  val namespace = "threestep"

  val inlineJScriptRequirements = Option(Array(
    Coproduct[Requirement](ShellCommandRequirement()),
      Coproduct[Requirement](InlineJavascriptRequirement())))

  val psOutputBinding = CommandOutputBinding(glob = Option(Coproduct[Glob]("ps-stdOut.txt")))

  val psOutputParameter =
    CommandOutputParameter(
      id = "ps-stdOut",
      `type` = Option(Coproduct(CwlType.File)),
      outputBinding = Option(psOutputBinding))

  val psClt = CommandLineTool(
    `class` = "CommandLineTool".narrow,
    outputs = Coproduct[CommandLineTool.Outputs](Array(psOutputParameter)),
    baseCommand = Option(Coproduct[BaseCommand]("ps")),
    stdout = Option(Coproduct[StringOrExpression]("ps-stdOut.txt")))

  val psWfStep  = WorkflowStep(
    id = Option("ps"),
    run = Coproduct[Run](psClt),
    out = Coproduct(Array("ps-stdOut")))

  val patternInput = CommandInputParameter(id = Option("pattern"), `type` = Option(Coproduct(CwlType.String)))

  val fileInput = CommandInputParameter(id = Option("file"), `type` = Option(Coproduct(CwlType.File)))

  def clb: String => CommandLineTool.Argument=
    s => Coproduct[Argument](CommandLineBinding(valueFrom = Option(Coproduct[StringOrExpression](s)), shellQuote  = Option(false)))

  val cgrepArgs = Option("grep $(inputs.pattern). $(inputs.file) | wc -l" split ' ' map clb)

  val cgrepOutputBinding = CommandOutputBinding(glob = Option(Coproduct[Glob]("cgrep-stdOut.txt")))

  val cgrepOutputParameter = CommandOutputParameter(id = "cgrep-stdOut", `type` = Option(Coproduct(CwlType.File)), outputBinding = Option(cgrepOutputBinding))

  val cgrepClt = CommandLineTool(
    inputs = Coproduct[CommandLineTool.Inputs](Array(patternInput, fileInput)),
    outputs = Coproduct(Array(cgrepOutputParameter)),
    `class` = "CommandLineTool".narrow,
    arguments = cgrepArgs,
    stdout = Option(Coproduct[StringOrExpression]("cgrep-stdOut.txt")),
    requirements = inlineJScriptRequirements)

  val patternCgrepWorkFlowStepInput = WorkflowStepInput(id = "pattern", source = Option(Coproduct("#pattern")))

  val fileCgrepWorkflowStepInput = WorkflowStepInput(id = "file", source = Option(Coproduct("ps/ps-stdOut")))

  val grepWfStep  = WorkflowStep(
    id = Option("cgrep"),
    in = Coproduct(Array(patternCgrepWorkFlowStepInput, fileCgrepWorkflowStepInput)),
    out = Coproduct[Outputs](Array(WorkflowStepOutput("cgrep-stdOut"))),
    run = Coproduct[Run](cgrepClt))

  val wcFileCommandInput = CommandInputParameter(
    id = Option("file"),
    `type` = Option(Coproduct(CwlType.File)))

  val wcArgs = Option("cat $(inputs.file) | wc -l" split ' ' map clb)

  val wcCltOutput = CommandOutputParameter(
    id = "wc-stdOut",
    `type` = Option(Coproduct(CwlType.File)),
    outputBinding = Option(CommandOutputBinding(glob = Option(Coproduct[Glob]("wc-stdOut.txt")))))

  val wcClt =
    CommandLineTool(
      `class` = "CommandLineTool".narrow,
      stdout = Option(Coproduct[StringOrExpression]("wc-stdOut.txt")),
      inputs = Coproduct(Array(wcFileCommandInput)),
      outputs = Coproduct(Array(wcCltOutput)),
      arguments = wcArgs,
      requirements = inlineJScriptRequirements)

  val wcWorkflowInput = WorkflowStepInput(
    id = "file",
    source = Option(Coproduct("ps/ps-stdOut")))

  val wcWorkflowStep = WorkflowStep(
    id = Option("wc"),
    in = Coproduct[WorkflowStep.Inputs](Array(wcWorkflowInput)),
    out = Coproduct[WorkflowStep.Outputs](Array(WorkflowStepOutput("wc-stdOut"))),
    run = Coproduct[WorkflowStep.Run](wcClt))

  val outputCgrep =
    WorkflowOutputParameter(
      id = Option("cgrep-stdOut"),
      `type` = Option(Coproduct[MyriadOutputType](CwlType.File)),
      outputSource = Option(Coproduct("#cgrep/cgrep-stdOut")))

  val outputWc =
    WorkflowOutputParameter(
      id = Option("wc-stdOut"),
      `type` = Option(Coproduct[MyriadOutputType](CwlType.File)),
      outputSource = Option(Coproduct("#wc/wc-stdOut")))

  val _outputs = Coproduct[WorkflowOutput](Array(outputCgrep, outputWc))

  val workflowPatternInput = InputParameter(id = Option("pattern"), `type` = Option(Coproduct[MyriadInputType](CwlType.String)))

  val _inputs = Coproduct[WorkflowInput](Array(workflowPatternInput))

  val threeStepWorkflow =
    new Workflow(
      inputs = _inputs,
      outputs = _outputs,
      steps = Coproduct[WorkflowSteps](Array(psWfStep, grepWfStep, wcWorkflowStep)))

  val yaml = encodeCwlWorkflow(threeStepWorkflow)

  println(yaml)
}
