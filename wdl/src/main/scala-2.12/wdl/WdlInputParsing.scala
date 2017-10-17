package wdl

import lenthall.Checked

import scala.util.Try
import wom.executable.Executable.InputParsingFunction
import wom.executable.Executable.ParsedInputMap
import wom.executable.Executable
import wom.types.WomType

private [wdl] object WdlInputParsing {
  private [wdl] lazy val inputCoercionFunction: InputParsingFunction = inputString => {
    import lenthall.validation.Checked._
    import lenthall.validation.Validation._
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
      womDefinition <- workflow.womDefinition.toEither
      executable <- Executable.withInputs(womDefinition, inputCoercionFunction, inputFile)
    } yield executable
  }
}
