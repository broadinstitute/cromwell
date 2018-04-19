package cwl

import cats.data.NonEmptyList
import shapeless.Witness
import wom.types.WomEnumerationType

trait EnumSchema {

  //This isn't seen in the spec but it has been observed in the conformance tests which makes one assume it is erroneously omitted
  val name: String
  val symbols: Array[String]
  val `type`: Witness.`"enum"`.T
  val label: Option[String]

  def toWomEnumerationType: WomEnumerationType = {
    val symbolIds = symbols.toList.map{
      s => s.substring(s.lastIndexOf("/") + 1)
    }
    WomEnumerationType(NonEmptyList.fromListUnsafe(symbolIds))
  }
}

case class InputEnumSchema(
  name: String,
  symbols: Array[String],
  `type`: Witness.`"enum"`.T = Witness("enum").value,
  label: Option[String] = None,
  inputBinding: Option[InputCommandLineBinding] = None) extends EnumSchema

case class OutputEnumSchema(
  name: String,
  symbols: Array[String],
  `type`: Witness.`"enum"`.T,
  label: Option[String],
  outputBinding: Option[CommandOutputBinding]) extends EnumSchema

