import wdl4s.values.WdlValue

import scala.util.Try

package object wdl4s {
  type WdlSource = String
  type WdlJson = String
  type WorkflowRawInputs = Map[FullyQualifiedName, Any]
  type WorkflowCoercedInputs = Map[FullyQualifiedName, WdlValue]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  type CallInputs = Map[String, WdlValue]
  type ImportResolver = String => WdlSource

  trait TsvSerializable {
    def tsvSerialize: Try[String]
  }

}
