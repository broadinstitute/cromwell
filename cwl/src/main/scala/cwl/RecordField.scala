package cwl

trait RecordField{
  val name: String
  val doc: Option[String]
}

case class OutputRecordField(
  name: String,
  `type`: MyriadOutputType,
  doc: Option[String],
  outputBinding: Option[CommandOutputBinding]) extends RecordField

case class InputRecordField(
  name: String,
  `type`: MyriadInputType,
  doc: Option[String],
  inputBinding: Option[InputCommandLineBinding],
  label: Option[String]) extends RecordField
