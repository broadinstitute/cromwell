package cromwell.binding.types

import cromwell.binding.WorkflowInput
import spray.json._

object WdlTypeJsonFormatter extends DefaultJsonProtocol {
  implicit object WdlTypeJsonFormat extends RootJsonFormat[WdlType] {
    def write(wdlType: WdlType) = JsString(wdlType.toWdlString)
    def read(value: JsValue) = ???
  }

  implicit object WorkflowInputJsonFormat extends RootJsonFormat[Seq[WorkflowInput]] {
    def write(input: Seq[WorkflowInput]) = {
      JsObject(input.map {input =>
        val optional = if (input.optional) "(optional) " else ""
        input.fqn -> JsString(s"$optional${input.types.map(_.toWdlString).mkString(", ")}")
      }.toMap)
    }
    def read(value: JsValue) = ???
  }
}

