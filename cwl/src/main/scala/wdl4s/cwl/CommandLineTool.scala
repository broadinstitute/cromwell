package wdl4s.cwl

import eu.timepit.refined._
import shapeless.{:+:, CNil}
import wdl4s.cwl.CwlType._

case class CommandInputParameter(
                                  id: String,
                                  label: Option[String] = None,
                                  secondaryFiles: Option[Array[ECMAScript :+: ECMAFunction :+: String :+: CNil]] = None,
                                  format: Option[ECMAScript :+: ECMAFunction :+: Array[String] :+: String :+: CNil] = None, //only valid when type: File
                                  streamable: Option[Boolean] = None, //only valid when type: File
                                  doc: Option[String :+: Array[String] :+: CNil] = None,
                                  inputBinding: Option[CommandLineBinding] = None,
                                  default: Option[CwlAny] = None,
                                  `type`: Option[MyriadInputType] = None)

case class CommandInputRecordSchema(
                                     `type`: W.`"record"`.T,
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
                                   `type`: W.`"enum"`.T,
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
                                    `type`: W.`"array"`.T,
                                    label: Option[String],
                                    inputBinding: Option[CommandLineBinding])


case class CommandOutputParameter(
                                   id: String,
                                   label: Option[String] = None,
                                   secondaryFiles: Option[ECMAScript :+: ECMAFunction :+: String :+: Array[ECMAScript :+: ECMAFunction :+: String :+: CNil] :+: CNil] = None,
                                   format: Option[ECMAScript :+: ECMAFunction :+: Array[String] :+: String :+: CNil] = None, //only valid when type: File
                                   streamable: Option[Boolean] = None, //only valid when type: File
                                   doc: Option[String :+: Array[String] :+: CNil] = None,
                                   outputBinding: Option[CommandOutputBinding] = None,
                                   `type`: Option[MyriadOutputType] = None)


