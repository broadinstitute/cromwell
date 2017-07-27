package wdl4s.cwl

import io.circe._
import io.circe.syntax._
import io.circe.yaml._
import io.circe.yaml.syntax._
import io.circe.generic.auto._
import io.circe.yaml.{parser => YamlParser}
import io.circe.parser._
import io.circe.shapes._
import io.circe.generic.auto._
import cats.syntax.either._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.string._
import eu.timepit.refined._
import io.circe.refined._
import io.circe.literal._

object CwlCodecs {

  def decodeCwl: Yaml => Either[Error, Cwl] = {
    import wdl4s.cwl.Implicits._

    YamlParser.
      parse(_).
      map(_.noSpaces).
      flatMap{json =>
        decode[CommandLineTool](json) orElse
          decode[Workflow](json)
      }
  }

  def encodeCwlCommandLineTool(commandLineTool: CommandLineTool): Json = {
    import io.circe.syntax._
    import wdl4s.cwl.Implicits.enumerationEncoder
    commandLineTool.asJson
  }

  def encodeCwlWorkflow(workflow: Workflow): Json = {
    import io.circe.syntax._
    import wdl4s.cwl.Implicits.enumerationEncoder
    workflow.asJson
  }

  def encodeCwl(cwl: Cwl): Json = {
    import io.circe.syntax._
    import wdl4s.cwl.Implicits.enumerationEncoder
    cwl match {
      case commandLineTool: CommandLineTool => commandLineTool.asJson
      case workflow: Workflow => workflow.asJson
    }
  }

  val jsonPrettyPrinter = io.circe.Printer.spaces2.copy(dropNullKeys = true, preserveOrder = true)
  val yamlPrettyPrinter = io.circe.yaml.Printer.spaces2.copy(dropNullKeys = true, preserveOrder = true)

  def cwlToJson(cwl: Cwl): String = jsonPrettyPrinter.pretty(encodeCwl(cwl))

  def cwlToYaml(cwl: Cwl): Yaml = yamlPrettyPrinter.pretty(encodeCwl(cwl))

}
