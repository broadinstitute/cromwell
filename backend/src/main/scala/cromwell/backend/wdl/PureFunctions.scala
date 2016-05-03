package cromwell.backend.wdl

import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.types.WdlStringType
import wdl4s.values.{WdlString, WdlValue}

import scala.util.{Failure, Success, Try}

trait PureFunctions { this: WdlStandardLibraryFunctions =>

  override def sub(params: Seq[Try[WdlValue]]): Try[WdlString] = {
    def extractArguments = params.size match {
      case 3 => Success((params.head, params(1), params(2)))
      case n => Failure(new IllegalArgumentException(s"Invalid number of parameters for engine function sub: $n. sub takes exactly 3 parameters."))
    }

    def validateArguments(values: (Try[WdlValue], Try[WdlValue], Try[WdlValue])) = values match {
      case (Success(strValue), Success(WdlString(pattern)), Success(replaceValue))
        if WdlStringType.isCoerceableFrom(strValue.wdlType) &&
          WdlStringType.isCoerceableFrom(replaceValue.wdlType) =>
        Success((strValue.valueString, pattern, replaceValue.valueString))
      case _ => Failure(new IllegalArgumentException(s"Invalid parameters for engine function sub: $values."))
    }

    for {
      args <- extractArguments
      (str, pattern, replace) <- validateArguments(args)
    } yield WdlString(pattern.r.replaceAllIn(str, replace))
  }
}
