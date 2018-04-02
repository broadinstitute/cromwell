package cwl
import shapeless.{Witness => W}

trait RecordSchema[T <: RecordField] {

  val name: String

  val `type`: W.`"record"`.T

  val fields: Option[Array[RecordField]]

  val label: Option[String]
}

case class InputRecordSchema(
  //NOTE Well: this does not actually appear in v1.0, but that is a  mistake and it should have been
  name: String,
  `type`: W.`"record"`.T,
  fields: Option[Array[InputRecordField]],
  label: Option[String]) extends RecordSchema[InputRecordField]

case class OutputRecordSchema(
  `type`: W.`"record"`.T,
  fields: Option[Array[OutputRecordField]],
  label: Option[String]) extends RecordSchema[OutputRecordField]
