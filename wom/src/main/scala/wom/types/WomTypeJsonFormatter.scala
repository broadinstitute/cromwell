package wom.types

import spray.json._
import wom.callable.Callable.InputDefinition
import wom.core.FullyQualifiedName

object WomTypeJsonFormatter extends DefaultJsonProtocol {
  implicit object WomTypeJsonFormat extends RootJsonFormat[WomType] {
    def write(womType: WomType) = JsString(womType.toDisplayString)
    def read(value: JsValue) = ???
  }

  implicit object WorkflowInputJsonFormat extends RootJsonFormat[Map[FullyQualifiedName, InputDefinition]] {
    def write(inputs: Map[FullyQualifiedName, InputDefinition]) = {
      JsObject(inputs map { case (fqn, input) =>
        val optional = if (input.optional) "(optional) " else ""
        fqn -> JsString(s"$optional${input.womType.toDisplayString}")
      })
    }
    def read(value: JsValue) = ???
  }
}

