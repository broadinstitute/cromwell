package cromwell.binding.types

import cromwell.binding.{FullyQualifiedName, WorkflowInput}
import spray.json._

object WdlTypeJsonFormatter extends DefaultJsonProtocol {
  implicit object WdlTypeJsonFormat extends RootJsonFormat[WdlType] {
    def write(wdlType: WdlType) = JsString(wdlType.toWdlString)
    def read(value: JsValue) = ???
  }

  implicit object WorkflowInputJsonFormat extends RootJsonFormat[Map[FullyQualifiedName, WorkflowInput]] {
    def write(inputs: Map[FullyQualifiedName, WorkflowInput]) = {
      JsObject(inputs map { case (fqn, input) =>
        val optional = if (input.optional) "(optional) " else ""
        fqn -> JsString(s"$optional${input.wdlType.toWdlString}")
      })
    }
    def read(value: JsValue) = ???
  }
}

