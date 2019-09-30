package cwl.encoder

import cwl._
import io.circe._

object CwlEncoder {

  def encodeCwl(cwl: Cwl): Json = cwl.fold(CwlEncoderFold)

  val jsonPrettyPrinter = io.circe.Printer.spaces2.copy(dropNullValues = true)

  def cwlToJson(cwl: Cwl): String = jsonPrettyPrinter.print(encodeCwl(cwl))
}


