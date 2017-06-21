package broad.cwl.model

import shapeless.{:+:, CNil}
import eu.timepit.refined.string._
import eu.timepit.refined._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.auto._
import io.circe.refined._
import enumeratum.Circe._

import enumeratum.Circe._
import io.circe.syntax._
import io.circe._
import io.circe.parser._
import io.circe.shapes._
import io.circe.generic.auto._
import shapeless.{:+:, CNil}
import cats._, implicits._//, instances._

import CWLVersion._

import io.circe.yaml.parser
import io.circe._

import eu.timepit.refined.string._
import eu.timepit.refined._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.auto._
import io.circe.refined._

case class CommandLineTool(
  inputs: CommandInputParameter :+: Map[CommandInputParameter#Id, CommandInputParameter#`type`] :+: Map[CommandInputParameter#Id, CommandInputParameter] :+: CNil,
  outputs: Array[CommandOutputParameter] :+: Map[CommandOutputParameter#Id, CommandOutputParameter#`type`] :+: Map[CommandOutputParameter#Id, CommandOutputParameter] :+: CNil,
  `class`: String,  //TODO: must be "CommandLineTool" or ?
  id: Option[String],
  requirements: Option[Array[Requirement]],
  hints: Option[Array[String]], //TODO: Any?
  label: Option[String],
  doc: Option[String],
  cwlVersion: Option[CWLVersion],
  baseCommand: Option[String :+: Array[String] :+: CNil],
  arguments: Option[Array[Expression :+: CommandLineBinding :+: String :+: CNil]],
  stdin: Option[Expression :+: String :+: CNil],
  stderr: Option[Expression :+: String :+: CNil],
  stdout: Option[Expression :+: String :+: CNil],
  successCodes: Option[Array[Int]],
  temporaryFailCodes: Option[Array[Int]],
  permanentFailCodes: Option[Array[Int]])

case class CommandInputParameter(
  id: Option[String],
  label: Option[String],
  secondaryFiles: Option[Array[Expression :+: String :+: CNil]],
  format: Option[Expression :+: Array[String] :+: String :+: CNil], //only valid when type: File
  streamable: Option[Boolean],//only valid when type: File
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
    CWLType :+:
    CommandInputRecordSchema :+:
    CommandInputEnumSchema :+:
    CommandInputArraySchema :+:
    String :+:
    Array[
      CWLType :+:
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
  secondaryFiles: Option[Expression :+: String :+: Array[Expression :+: String :+: CNil] :+: CNil],
  format: Expression :+: Array[String] :+: String :+: CNil, //only valid when type: File
  streamable: Option[Boolean],//only valid when type: File
  doc: Option[String :+: Array[String] :+: CNil],
  outputBinding: Option[CommandOutputBinding],
  `type`: MyriadOutputType) {

  type `type` = String
  type Id = String
}

