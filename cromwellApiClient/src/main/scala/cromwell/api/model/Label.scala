package cromwell.api.model

import spray.json.{DefaultJsonProtocol, JsObject, JsString, JsValue, RootJsonFormat}

object LabelsJsonFormatter extends DefaultJsonProtocol {
  implicit object LabelJsonFormat extends RootJsonFormat[List[Label]] {
    def write(l: List[Label]) = JsObject(l map { label => label.key -> JsString(label.value)} :_* )
    def read(value: JsValue) = value match {
      case JsObject(x) => x map { case (k, JsString(v)) => Label(k, v) } toList
    }
  }
}

case class Label(key: String, value: String)
