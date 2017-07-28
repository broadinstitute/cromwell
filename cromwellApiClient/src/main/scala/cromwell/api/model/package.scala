package cromwell.api

import java.time.OffsetDateTime

import spray.json.{DefaultJsonProtocol, JsString, JsValue, RootJsonFormat}

package object model {

  implicit val OffsetDateTimeJsonFormat = OffsetDateTimeJsonFormatter.OffsetDateTimeFormat

  object OffsetDateTimeJsonFormatter extends DefaultJsonProtocol {
    object OffsetDateTimeFormat extends RootJsonFormat[OffsetDateTime] {
      def write(odt: OffsetDateTime) = new JsString(odt.toString)
      def read(value: JsValue) = value match {
        case JsString(string) => OffsetDateTime.parse(string)
        case other => throw new UnsupportedOperationException(s"Cannot deserialize $other into an OffsetDateTime")
      }
    }
  }
}
