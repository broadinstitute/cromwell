package wes2cromwell

import spray.json.{ DefaultJsonProtocol, JsObject, JsonFormat, JsonParser }

case class CromwellCallsMetadata(
  shardIndex: Option[Int],
  returnCode: Option[Int],
  start: Option[String],
  end: Option[String],
  stdout: Option[String],
  stderr: Option[String]
)

object CromwellCallsMetadata {
  import DefaultJsonProtocol._
  implicit val cromwellCallsMetadataFormat: JsonFormat[CromwellCallsMetadata] = jsonFormat6(CromwellCallsMetadata.apply)
}

case class CromwellSubmittedFiles(
  workflow: Option[String],
  workflowType: String,
  workflowTypeVersion: String,
  options: Option[String],
  inputs: Option[String],
  labels: Option[String]
)

object CromwellSubmittedFiles {
  import DefaultJsonProtocol._
  implicit val cromwellSubmittedFilesFormat: JsonFormat[CromwellSubmittedFiles] = jsonFormat6(CromwellSubmittedFiles.apply)
}

case class CromwellMetadata(
  workflowName: Option[String],
  id: String,
  status: String,
  start: Option[String],
  end: Option[String],
  submittedFiles: CromwellSubmittedFiles,
  outputs: Option[JsObject],
  calls: Option[Map[String, Seq[CromwellCallsMetadata]]]
)

object CromwellMetadata {
  import DefaultJsonProtocol._
  implicit val cromwellMetadataFormat: JsonFormat[CromwellMetadata] = jsonFormat8(CromwellMetadata.apply)

  def toCromwellMetadata(json: String): CromwellMetadata = {
    val jsonAst = JsonParser(json)
    jsonAst.convertTo[CromwellMetadata]
  }
}
