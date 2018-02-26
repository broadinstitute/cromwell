package wdl.shared.transforms.wdlom2wom

// TODO 2.11 remove this "import cats.syntax.either._"
import common.Checked
import common.validation.Checked._
import common.collections.EnhancedCollections._
import wom.callable.WorkflowDefinition
import wom.core.WorkflowJson
import wom.executable.{Executable, WomBundle}
import wom.executable.Executable.{InputParsingFunction, ParsedInputMap}
import wom.types.WomType

import scala.util.Try

object WdlSharedInputParsing {
  private lazy val inputCoercionFunction: InputParsingFunction = inputString => {
    import common.validation.Checked._
    import common.validation.Validation._
    import spray.json._

    Try(inputString.parseJson).toErrorOr.toEither flatMap {
      case JsObject(fields) => fields.map({
        case (key, jsValue) => key -> { womType: WomType => womType.coerceRawValue(jsValue).toErrorOr }
      }).validNelCheck
      case other => s"WDL input file must be a valid Json object. Found a ${other.getClass.getSimpleName}".invalidNelCheck[ParsedInputMap]
    }
  }

  def buildWomExecutable(bundle: WomBundle, inputs: Option[WorkflowJson]): Checked[Executable] = {

    val workflows: Set[WorkflowDefinition] = bundle.callables.filterByType[WorkflowDefinition]
    val executableWorkflowCheck = if (workflows.size == 1) {
      workflows.head.validNelCheck
    } else {
      s"Cannot convert WOM bundle to executable: required exactly one workflow definition but got ${workflows.size}".invalidNelCheck
    }

    for {
      executableWorkflow <- executableWorkflowCheck
      executable <- Executable.withInputs(executableWorkflow, inputCoercionFunction, inputs)
    } yield executable
  }
}
