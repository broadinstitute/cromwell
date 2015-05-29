package cromwell.binding.types

import spray.json._

object WdlTypeJsonFormatter extends DefaultJsonProtocol {
  implicit object WdlTypeJsonFormat extends RootJsonFormat[WdlType] {
    def write(wdlType: WdlType) = JsString(wdlType.toWdlString)
    def read(value: JsValue) = ???
  }
}

