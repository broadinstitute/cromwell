package cromwell.backend.wdl

import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.types.{WdlArrayType, WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlArray, WdlFile, WdlFloat, WdlInteger, WdlString, WdlValue}

import scala.util.{Failure, Success, Try}

case object OnlyPureFunctions extends WdlStandardLibraryFunctions with PureFunctions {
  override def readFile(path: String): String = throw new NotImplementedError("readFile not available in PureNoFunctions.")
  override def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue] = throw new NotImplementedError("read_json not available in PureNoFunctions.")
  override def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError("write_json not available in PureNoFunctions.")
  override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = throw new NotImplementedError("size not available in PureNoFunctions.")
  override def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError("write_tsv not available in PureNoFunctions.")
  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError("stdout not available in PureNoFunctions.")
  override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedError("glob not available in PureNoFunctions.")
  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = throw new NotImplementedError("writeTempFile not available in PureNoFunctions.")
  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError("stderr not available in PureNoFunctions.")
}

trait PureFunctions { this: WdlStandardLibraryFunctions =>

  def range(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    def extractAndValidateArguments = params.size match {
      case 1 => validateArguments(params.head)
      case n => Failure(new IllegalArgumentException(s"Invalid number of parameters for engine function seq: $n. Ensure seq(x: WdlInteger) takes exactly 1 parameters."))
    }

    def validateArguments(value: Try[WdlValue]) = value match {
      case Success(intValue: WdlValue) if WdlIntegerType.isCoerceableFrom(intValue.wdlType) =>
        Integer.valueOf(intValue.valueString) match {
          case i if i >= 0 => Success(i)
          case n => Failure(new IllegalArgumentException(s"Parameter to seq must be greater than or equal to 0 (but got $n)"))
        }
      case _ => Failure(new IllegalArgumentException(s"Invalid parameter for engine function seq: $value."))
    }

    extractAndValidateArguments map { intValue => WdlArray(WdlArrayType(WdlIntegerType), (0 until intValue).map(WdlInteger(_))) }
  }

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
