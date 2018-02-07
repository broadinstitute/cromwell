package wdl.transforms.draft2.wdlom2wom

// TODO 2.11 remove this "import cats.syntax.either._"
import cats.syntax.either._
import common.Checked
import wdl.WdlWorkflow

import scala.util.Try
import wom.executable.Executable.InputParsingFunction
import wom.executable.Executable.ParsedInputMap
import wom.executable.Executable
import wom.types.WomType
import wom.transforms.WomWorkflowDefinitionMaker.ops._

object WdlDraft2InputParsing {
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

  def buildWomExecutable(workflow: WdlWorkflow, inputFile: Option[String]): Checked[Executable] = {
    for {
      womDefinition <- workflow.toWomWorkflowDefinition.toEither
      executable <- Executable.withInputs(womDefinition, inputCoercionFunction, inputFile)
    } yield executable
  }
}
