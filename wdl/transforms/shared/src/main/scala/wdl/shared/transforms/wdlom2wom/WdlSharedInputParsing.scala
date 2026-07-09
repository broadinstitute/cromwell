package wdl.shared.transforms.wdlom2wom

import common.Checked
import wom.core.WorkflowJson
import wom.executable.Executable.{InputParsingFunction, ParsedInputMap}
import wom.executable.{Executable, WomBundle}
import wom.expression.IoFunctionSet
import wom.types.WomType

import scala.util.Try

object WdlSharedInputParsing {

  // "my_wf.my_task.runtime.cpu" -> ("my_wf.my_task.runtime", "cpu")
  final val runtimeOverrideRegex = "(.*\\.runtime)\\.([^\\.]*)$".r

  private lazy val inputCoercionFunction: InputParsingFunction = inputString => {
    import common.validation.Checked._
    import common.validation.Validation._
    import spray.json._

    Try(inputString.parseJson).toErrorOr.toEither flatMap {
      case JsObject(fields) =>
        // Accumulate runtime overrides, we want to output them as a single object rather than individual fields
        // For example, this input JSON string:
        // {
        //   "task1.runtime.cpu": "2",
        //   "task1.runtime.memory": "4 GB",
        //   "task1.input1": "some value",
        //   "task2.runtime.disk": "10 GB"
        // }
        // Should result in a ParsedInputMap with this structure:
        // {
        //   "task1.runtime": {
        //     "cpu": "2",
        //     "memory": "4 GB"
        //   },
        //   "task1.input1": "some value",
        //   "task2.runtime": {
        //     "disk": "10 GB"
        //   }
        // }
        var runtimeOverrideFields: Map[String, Map[String, JsValue]] = Map.empty

        val inputFields = fields.flatMap { case (key, jsValue) =>
          runtimeOverrideRegex.findFirstMatchIn(key) match {
            case Some(m) =>
              val runtimeOverrideSetKey = m.group(1) // e.g. "task1.runtime"
              val runtimeOverrideAttributeKey = m.group(2) // e.g. "cpu"
              val existingAttributes = runtimeOverrideFields.getOrElse(runtimeOverrideSetKey, Map.empty)
              val updatedAttributes = existingAttributes + (runtimeOverrideAttributeKey -> jsValue)
              runtimeOverrideFields = runtimeOverrideFields + (runtimeOverrideSetKey -> updatedAttributes)
              None
            case None =>
              Option(key -> { womType: WomType => womType.coerceRawValue(jsValue).toErrorOr })
          }
        }

        val runtimeOverrides = runtimeOverrideFields.map { case (runtimeKey, attributesMap) =>
          val jsObject = JsObject(attributesMap)
          runtimeKey -> { womType: WomType =>
            womType.coerceRawValue(jsObject).toErrorOr
          }
        }

        (inputFields ++ runtimeOverrides).validNelCheck
      case other =>
        s"WDL input file must be a valid Json object. Found a ${other.getClass.getSimpleName}"
          .invalidNelCheck[ParsedInputMap]
    }
  }

  // Used for unit testing
  def parseInputs(inputString: String): Checked[ParsedInputMap] =
    inputCoercionFunction(inputString)

  def buildWomExecutable(bundle: WomBundle,
                         inputs: Option[WorkflowJson],
                         ioFunctions: IoFunctionSet,
                         strictValidation: Boolean
  ): Checked[Executable] =
    for {
      ec <- bundle.toExecutableCallable
      executable <- Executable.withInputs(ec, inputCoercionFunction, inputs, ioFunctions, strictValidation)
    } yield executable
}
