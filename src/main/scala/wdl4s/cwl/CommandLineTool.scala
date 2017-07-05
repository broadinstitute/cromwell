package wdl4s.cwl

import eu.timepit.refined._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.string._
import shapeless.{:+:, CNil}
import wdl4s.cwl.CwlVersion._
import wdl4s.cwl.CwlType._


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
                                        CNil] :+:
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

