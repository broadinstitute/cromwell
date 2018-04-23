package cwl

import cats.data.NonEmptyList
import shapeless.Witness
import wom.types.WomEnumerationType

trait EnumSchema {

  val name: Option[String]
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
  name: Option[String],
  symbols: Array[String],
  `type`: Witness.`"enum"`.T = Witness("enum").value,
  label: Option[String] = None,
  inputBinding: Option[InputCommandLineBinding] = None) extends EnumSchema

case class OutputEnumSchema(
  name: Option[String],
  symbols: Array[String],
  `type`: Witness.`"enum"`.T,
  label: Option[String],
  outputBinding: Option[CommandOutputBinding]) extends EnumSchema

