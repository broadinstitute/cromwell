package cwl

import cats.data.NonEmptyList
import io.circe._
import io.circe.parser._
import io.circe.shapes._
import io.circe.generic.auto._
import eu.timepit.refined.string._
import io.circe.refined._
import io.circe.literal._
import shapeless.Coproduct
import cats.syntax.either._
import cats.syntax.show._


object CwlCodecs {
  import Implicits._

  def decodeCwl(in: String): Either[NonEmptyList[String], Cwl] = {
    //try to parse both and combine errors if they fail
    (decode[Workflow](in), decode[CommandLineTool](in)) match {
      case (Right(wf), _) => Coproduct[Cwl](wf).asRight
      case (_, Right(clt)) => Coproduct[Cwl](clt).asRight
      case (Left(wfError), Left(cltError)) =>
        //This is not really suppressed but there is no other way to compose errors at the Exception level API AFAIK
        NonEmptyList.of(
          s"Workflow parsing error: ${wfError.show}",
          s"Command Line Tool parsing error: ${cltError.show}"
        ).asLeft
    }
  }

  def encodeCwlCommandLineTool(commandLineTool: CommandLineTool): Json = {
    import io.circe.syntax._
    import cwl.Implicits.enumerationEncoder
    commandLineTool.asJson
  }

  def encodeCwlWorkflow(workflow: Workflow): Json = {
    import io.circe.syntax._
    import cwl.Implicits.enumerationEncoder
    workflow.asJson
  }

  def encodeCwl(cwl: Cwl): Json = cwl.fold(CwlEncoder)

  val jsonPrettyPrinter = io.circe.Printer.spaces2.copy(dropNullValues = true, preserveOrder = true)

  def cwlToJson(cwl: Cwl): String = jsonPrettyPrinter.pretty(encodeCwl(cwl))
}
