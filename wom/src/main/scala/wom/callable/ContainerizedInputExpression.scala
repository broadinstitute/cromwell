package wom.callable

import common.validation.ErrorOr.ErrorOr
import wom.expression.IoFunctionSet
import wom.values.WomValue

/**
  * This is an expression that uses "containerized" (a.k.a. "mapped" input values) in its evaluation.
  */
trait ContainerizedInputExpression {
  def evaluate(hostInputValues: Map[String, WomValue],
               containerizedInputValues: Map[String, WomValue],
               ioFunctionSet: IoFunctionSet): ErrorOr[AdHocValue]
}

/**
  * Tracks womValues that should optionally update an input from a Yaml/Json.
  *
  * Example:
  * https://github.com/common-workflow-language/common-workflow-language/commit/a486e91f966780e6d2321ecaef57f3d610f9f632#diff-cf75ec4d6999e932fb716212a41c8f73R6
  *
  * The input Yaml v1.0/wc-job.json says that the input is `whale.txt`.
  * However the CWL v1.0/initialwork-path.cwl says that the file should be changed to `bob.txt`.
  * This change happens during the evaluation of the expression `$(inputs.file1)`.
  */
final case class AdHocValue(womValue: WomValue, mutableInputOption: Option[String])
