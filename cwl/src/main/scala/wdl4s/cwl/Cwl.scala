package wdl4s.cwl

import shapeless.{:+:, CNil, Coproduct, Witness}
import shapeless.syntax.singleton._
import eu.timepit.refined._
import CwlVersion._
import wdl4s.cwl.CommandLineTool.{apply => _, _}


sealed trait Cwl {

  val cwlVersion: Option[CwlVersion]
}

case class Workflow(
  cwlVersion: Option[CwlVersion] = Option(CwlVersion.Version1),
  `class` : Workflow.`class`.type = Workflow.`class`,
  inputs: Array[InputParameter] = Array.empty,
  outputs: Array[WorkflowOutputParameter] = Array.empty,
  steps: Array[WorkflowStep]) extends Cwl

object Workflow {
  val `class` : Witness.`"Workflow"`.T = "Workflow".narrow
}

/**
  *
  * @param inputs
  * @param outputs
  * @param `class` This _should_ always be "CommandLineTool," however the spec does not -er- specify this.
  * @param id
  * @param requirements
  * @param hints
  * @param label
  * @param doc
  * @param cwlVersion
  * @param baseCommand
  * @param arguments
  * @param stdin
  * @param stderr
  * @param stdout
  * @param successCodes
  * @param temporaryFailCodes
  * @param permanentFailCodes
  */
case class CommandLineTool(
                            inputs: Array[CommandInputParameter] = Array.empty,
                            outputs: Array[CommandOutputParameter] = Array.empty,
                            `class`: CommandLineTool.`class`.type = CommandLineTool.`class`,
                            id: Option[String] = None,
                            requirements: Option[Array[Requirement]] = None,
                            hints: Option[Array[String]] = None, //TODO: Any?
                            label: Option[String] = None,
                            doc: Option[String] = None,
                            cwlVersion: Option[CwlVersion] = Option(CwlVersion.Version1),
                            baseCommand: Option[BaseCommand] = None,
                            arguments: Option[Array[Argument]] = None,
                            stdin: Option[StringOrExpression] = None,
                            stderr: Option[StringOrExpression] = None,
                            stdout: Option[StringOrExpression] = None,
                            successCodes: Option[Array[Int]] = None,
                            temporaryFailCodes: Option[Array[Int]] = None,
                            permanentFailCodes: Option[Array[Int]] = None) extends Cwl

object CommandLineTool {
  val `class` : Witness.`"CommandLineTool"`.T = "CommandLineTool".narrow

  type StringOrExpression = ECMAScriptExpression :+: String :+: CNil

  type BaseCommand = String :+: Array[String] :+: CNil

  type Argument = ECMAScriptExpression :+: CommandLineBinding :+: String :+: CNil
}
