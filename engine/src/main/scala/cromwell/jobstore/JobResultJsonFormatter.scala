package cromwell.jobstore

import cromwell.util.JsonFormatting.WomValueJsonFormatter.WomValueJsonFormat
import cromwell.core.CallOutputs
import spray.json.{DefaultJsonProtocol, JsValue, RootJsonFormat, _}
import wom.JobOutput

object JobResultJsonFormatter extends DefaultJsonProtocol {
  implicit object ThrowableFormat extends RootJsonFormat[Throwable] {
    def write(value: Throwable) = value.getMessage.toJson
    def read(value: JsValue): Throwable = new Exception(value.convertTo[String])
  }

  implicit object JobOutputFormat extends RootJsonFormat[JobOutput] {
    def write(value: JobOutput) = value.womValue.toJson
    def read(value: JsValue): JobOutput = JobOutput(WomValueJsonFormat.read(value))
  }

  implicit object CallOutputsFormat extends RootJsonFormat[CallOutputs] {
    def write(value: CallOutputs) = value.outputs.map({case (port, v) => port.identifier.fullyQualifiedName.value -> v }).toJson
    def read(value: JsValue): CallOutputs = throw new UnsupportedOperationException("Cannot deserialize outputs to output ports")
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
