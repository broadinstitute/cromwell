package wdl4s.cwl

import shapeless.{:+:, CNil}
import eu.timepit.refined._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.string._
import CwlVersion._
import CwlType._

sealed trait Cwl {

  val cwlVersion: Option[CwlVersion]

  val `class`: String
}

case class Workflow(
  cwlVersion: Option[CwlVersion] = None,
  `class`: String,
  inputs: WorkflowInput,
  outputs: WorkflowOutput,
  steps: WorkflowSteps) extends Cwl

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
                            inputs:
                              CommandInputParameter :+:
                              Map[CommandInputParameter#Id, CommandInputParameter#`type`] :+:
                              Map[CommandInputParameter#Id, CommandInputParameter] :+:
                              CNil,
                            outputs:
                              Array[CommandOutputParameter] :+:
                              Map[CommandOutputParameter#Id, CommandOutputParameter#`type`] :+:
                              Map[CommandOutputParameter#Id, CommandOutputParameter] :+:
                              CNil,
                            `class`: String,
                            id: Option[String],
                            requirements: Option[Array[Requirement]],
                            hints: Option[Array[String]], //TODO: Any?
                            label: Option[String],
                            doc: Option[String],
                            cwlVersion: Option[CwlVersion],
                            baseCommand: Option[String :+: Array[String] :+: CNil],
                            arguments: Option[Array[ECMAScriptExpression :+: CommandLineBinding :+: String :+: CNil]],
                            stdin: Option[ECMAScriptExpression :+: String :+: CNil],
                            stderr: Option[ECMAScriptExpression :+: String :+: CNil],
                            stdout: Option[ECMAScriptExpression :+: String :+: CNil],
                            successCodes: Option[Array[Int]],
                            temporaryFailCodes: Option[Array[Int]],
                            permanentFailCodes: Option[Array[Int]]) extends Cwl
