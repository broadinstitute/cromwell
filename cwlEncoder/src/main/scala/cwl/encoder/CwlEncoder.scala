package cwl.encoder

import cwl._
import io.circe._

object CwlEncoder {

  def encodeCwl(cwl: Cwl): Json = cwl.fold(CwlEncoderFold)

  val jsonPrettyPrinter = io.circe.Printer.spaces2.copy(dropNullValues = true, preserveOrder = true)

  def cwlToJson(cwl: Cwl): String = jsonPrettyPrinter.pretty(encodeCwl(cwl))
}


