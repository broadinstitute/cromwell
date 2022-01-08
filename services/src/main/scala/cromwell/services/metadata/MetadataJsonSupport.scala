package cromwell.services.metadata

import java.time.OffsetDateTime

import cromwell.core.WorkflowId
import spray.json.{DefaultJsonProtocol, JsString, JsValue, RootJsonFormat}
import common.util.TimeUtil._

object MetadataJsonSupport extends DefaultJsonProtocol {
  implicit object WorkflowIdJsonFormatter extends RootJsonFormat[WorkflowId] {
    def write(id: WorkflowId) = JsString(id.id.toString)
    def read(value: JsValue) = throw new UnsupportedOperationException("Reading WorkflowId from JSON is currently unsupported")
  }

  implicit object MetadataTypeFormatter extends RootJsonFormat[MetadataType] {
    def write(t: MetadataType) = JsString(t.typeName)
    def read(value: JsValue) = throw new UnsupportedOperationException("Reading MetadataType from JSON is currently unsupported")
  }

  implicit object OffsetDateTimeFormatter extends RootJsonFormat[OffsetDateTime] {
    def write(offsetDateTime: OffsetDateTime) = new JsString(offsetDateTime.toUtcMilliString)
    def read(value: JsValue) = throw new UnsupportedOperationException("Reading OffsetDateTime from JSON is currently unsupported")
  }

  implicit val MetadataValueFormat = jsonFormat2(MetadataValue.apply)
  implicit val MetadataJobKeyFormat = jsonFormat3(MetadataJobKey)
  implicit val MetadataKeyFormat = jsonFormat3(MetadataKey.apply)
  implicit val MetadataEventFormat = jsonFormat3(MetadataEvent.apply)
}
