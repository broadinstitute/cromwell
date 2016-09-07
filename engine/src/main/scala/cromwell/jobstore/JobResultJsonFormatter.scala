package cromwell.jobstore

import cromwell.core.JobOutput
import cromwell.util.JsonFormatting.WdlValueJsonFormatter
import WdlValueJsonFormatter.WdlValueJsonFormat
import spray.json.{DefaultJsonProtocol, JsValue, RootJsonFormat, _}

object JobResultJsonFormatter extends DefaultJsonProtocol {
  implicit object ThrowableFormat extends RootJsonFormat[Throwable] {
    def write(value: Throwable) = value.getMessage.toJson
    def read(value: JsValue): Throwable = new Exception(value.convertTo[String])
  }

  implicit object JobOutputFormat extends RootJsonFormat[JobOutput] {
    def write(value: JobOutput) = value.wdlValue.toJson
    def read(value: JsValue): JobOutput = JobOutput(WdlValueJsonFormat.read(value))
  }

  implicit val JobResultSuccessFormat = jsonFormat2(JobResultSuccess)
  implicit val JobResultFailureFormat = jsonFormat3(JobResultFailure)

  implicit object JobResultFormat extends RootJsonFormat[JobResult] {
    def write(value: JobResult) = JsObject(
      value.getClass.getSimpleName -> (value match {
        case x: JobResultSuccess => x.toJson
        case x: JobResultFailure => x.toJson
      })
    )
    def read(value: JsValue): JobResult = {
      val fields = value.asJsObject.fields
      if (fields.contains("JobResultSuccess")) {
        fields("JobResultSuccess").convertTo[JobResultSuccess]
      } else {
        fields("JobResultFailure").convertTo[JobResultFailure]
      }
    }
  }
}