package cromwell.binding.values

import spray.json._

object WdlValueJsonFormatter extends DefaultJsonProtocol {
  implicit object WdlValueJsonFormat extends RootJsonFormat[WdlValue] {
    def write(value: WdlValue) = value match {
      case s:WdlString => JsString(s.value)
      case i:WdlInteger => JsNumber(i.value)
      case f:WdlFloat => JsNumber(f.value)
      case b:WdlBoolean => JsBoolean(b.value)
      case f:WdlFile => JsString(f.value.toAbsolutePath.toString)
      case o:WdlObject => JsObject()
    }
    def read(value: JsValue) = ???
  }
}
