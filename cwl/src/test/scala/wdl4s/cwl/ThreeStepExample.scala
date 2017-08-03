package wdl4s.cwl

import shapeless._
import syntax.singleton._
import wdl4s.cwl.CommandLineTool.{Argument, BaseCommand, StringOrExpression}
import wdl4s.cwl.CommandOutputBinding.Glob
import wdl4s.cwl.WorkflowStep.{Outputs, Run}
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
    outputs = Array(psOutputParameter),
    baseCommand = Option(Coproduct[BaseCommand]("ps")),
    stdout = Option(Coproduct[StringOrExpression]("ps-stdOut.txt")))

  val psWfStep  = WorkflowStep(
    id = "ps",
    run = Coproduct[Run](psClt),
    out = Coproduct(Array("ps-stdOut")))

  val patternInput = CommandInputParameter(id = "pattern", `type` = Option(Coproduct(CwlType.String)))

  val fileInput = CommandInputParameter(id = "file", `type` = Option(Coproduct(CwlType.File)))

  def clb: String => CommandLineTool.Argument=
    s => Coproduct[Argument](CommandLineBinding(valueFrom = Option(Coproduct[StringOrExpression](s)), shellQuote  = Option(false)))

  val cgrepArgs = Option("grep $(inputs.pattern). $(inputs.file) | wc -l" split ' ' map clb)

  val cgrepOutputBinding = CommandOutputBinding(glob = Option(Coproduct[Glob]("cgrep-stdOut.txt")))

  val cgrepOutputParameter = CommandOutputParameter(id = "cgrep-stdOut", `type` = Option(Coproduct(CwlType.File)), outputBinding = Option(cgrepOutputBinding))

  val cgrepClt = CommandLineTool(
    inputs = Array(patternInput, fileInput),
    outputs = Array(cgrepOutputParameter),
    `class` = "CommandLineTool".narrow,
    arguments = cgrepArgs,
    stdout = Option(Coproduct[StringOrExpression]("cgrep-stdOut.txt")),
    requirements = inlineJScriptRequirements)

  val patternCgrepWorkFlowStepInput = WorkflowStepInput(id = "pattern", source = Option(Coproduct("#pattern")))

  val fileCgrepWorkflowStepInput = WorkflowStepInput(id = "file", source = Option(Coproduct("ps/ps-stdOut")))

  val grepWfStep  = WorkflowStep(
    id = "cgrep",
    in = Array(patternCgrepWorkFlowStepInput, fileCgrepWorkflowStepInput),
    out = Coproduct[Outputs](Array(WorkflowStepOutput("cgrep-stdOut"))),
    run = Coproduct[Run](cgrepClt))

  val wcFileCommandInput = CommandInputParameter(
    id = "file",
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
      inputs = Array(wcFileCommandInput),
      outputs = Array(wcCltOutput),
      arguments = wcArgs,
      requirements = inlineJScriptRequirements)

  val wcWorkflowInput = WorkflowStepInput(
    id = "file",
    source = Option(Coproduct("ps/ps-stdOut")))

  val wcWorkflowStep = WorkflowStep(
    id = "wc",
    in = Array(wcWorkflowInput),
    out = Coproduct[WorkflowStep.Outputs](Array(WorkflowStepOutput("wc-stdOut"))),
    run = Coproduct[WorkflowStep.Run](wcClt))

  val outputCgrep =
    WorkflowOutputParameter(
      id = "cgrep-stdOut",
      `type` = Option(Coproduct[MyriadOutputType](CwlType.File)),
      outputSource = Option(Coproduct("#cgrep/cgrep-stdOut")))

  val outputWc =
    WorkflowOutputParameter(
      id = "wc-stdOut",
      `type` = Option(Coproduct[MyriadOutputType](CwlType.File)),
      outputSource = Option(Coproduct("#wc/wc-stdOut")))

  val _outputs = Array(outputCgrep, outputWc)

  val workflowPatternInput = InputParameter(id = "pattern", `type` = Option(Coproduct[MyriadInputType](CwlType.String)))

  val _inputs = Array(workflowPatternInput)

  val threeStepWorkflow =
    new Workflow(
      inputs = _inputs,
      outputs = _outputs,
      steps = Array(psWfStep, grepWfStep, wcWorkflowStep))

  val yaml = encodeCwlWorkflow(threeStepWorkflow)

  println(yaml)
}
