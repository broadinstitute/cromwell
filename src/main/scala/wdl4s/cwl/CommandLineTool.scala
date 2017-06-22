package wdl4s.cwl

import eu.timepit.refined._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.string._
import shapeless.{:+:, CNil}
import wdl4s.cwl.CwlVersion._
import wdl4s.cwl.CwlType._

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
                            permanentFailCodes: Option[Array[Int]])

case class CommandInputParameter(
                                  id: Option[String],
                                  label: Option[String],
                                  secondaryFiles: Option[Array[ECMAScriptExpression :+: String :+: CNil]],
                                  format: Option[ECMAScriptExpression :+: Array[String] :+: String :+: CNil], //only valid when type: File
                                  streamable: Option[Boolean], //only valid when type: File
                                  doc: Option[String :+: Array[String] :+: CNil],
                                  inputBinding: Option[CommandLineBinding],
                                  default: Option[String], //TODO Any type here
                                  `type`: Option[MyriadInputType]

                                ) {
  type Id = String
  type `type` = MyriadCommandInputType
}

case class CommandInputRecordSchema(
                                     `type`: String Refined MatchesRegex[W.`"record"`.T],
                                     fields: Option[Array[CommandInputRecordField]],
                                     label: Option[String])

case class CommandInputRecordField(
                                    name: String,
                                    `type`: MyriadInputType,
                                    doc: Option[String],
                                    inputBinding: Option[CommandLineBinding],
                                    label: Option[String])

case class CommandInputEnumSchema(
                                   symbols: Array[String],
                                   `type`: String Refined MatchesRegex[W.`"enum"`.T],
                                   label: Option[String],
                                   inputBinding: Option[CommandLineBinding])

case class CommandInputArraySchema(
                                    items:
                                    CwlType :+:
                                      CommandInputRecordSchema :+:
                                      CommandInputEnumSchema :+:
                                      CommandInputArraySchema :+:
                                      String :+:
                                      Array[
                                        CwlType :+:
                                          CommandInputRecordSchema :+:
                                          CommandInputEnumSchema :+:
                                          CommandInputArraySchema :+:
                                          String :+:
                                          CNil
                                        ] :+:
                                      CNil,
                                    `type`: String Refined MatchesRegex[W.`"array"`.T],
                                    label: Option[String],
                                    inputBinding: Option[CommandLineBinding])


case class CommandOutputParameter(
                                   id: String,
                                   label: Option[String],
                                   secondaryFiles: Option[ECMAScriptExpression :+: String :+: Array[ECMAScriptExpression :+: String :+: CNil] :+: CNil],
                                   format: ECMAScriptExpression :+: Array[String] :+: String :+: CNil, //only valid when type: File
                                   streamable: Option[Boolean], //only valid when type: File
                                   doc: Option[String :+: Array[String] :+: CNil],
                                   outputBinding: Option[CommandOutputBinding],
                                   `type`: MyriadOutputType) {

  type `type` = String
  type Id = String
}

