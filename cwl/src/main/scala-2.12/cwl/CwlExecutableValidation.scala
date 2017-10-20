package cwl

import io.circe._
import io.circe.Decoder
import io.circe.shapes._
import io.circe.generic.auto._
import io.circe.refined._
import io.circe.yaml
import io.circe.literal._
import common.Checked
import common.validation.Checked._
import wom.callable.Callable
import wom.executable.Executable
import wom.executable.Executable.{InputParsingFunction, ParsedInputMap}

// See explanation as to why there are 2 versions of this in ExecutableValidation
object CwlExecutableValidation {

  implicit val f = implicitly[Decoder[File]]

  // Decodes the input file, and build the ParsedInputMap
  private val inputCoercionFunction: InputParsingFunction =
    inputFile => {
      yaml.parser.parse(inputFile).flatMap(_.as[Map[String, MyriadInputValue]]) match {
        case Left(error) => error.getMessage.invalidNelCheck[ParsedInputMap]
        case Right(inputValue) => inputValue.map({ case (key, value) => key -> value.fold(CwlInputCoercion) }).validNelCheck
      }
    }

  def builWomExecutable(callable: Checked[Callable], inputFile: Option[String]): Checked[Executable] = {
    for {
      womDefinition <- callable
      executable <- Executable.withInputs(womDefinition, inputCoercionFunction, inputFile)
    } yield executable
  }
}
