package cromwell.backend.wdl

import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.types.{WdlArrayType, WdlIntegerType, WdlStringType, WdlType}
import wdl4s.values.{WdlArray, WdlFile, WdlFloat, WdlInteger, WdlString, WdlValue}

import scala.util.{Failure, Success, Try}

case object OnlyPureFunctions extends WdlStandardLibraryFunctions with PureFunctions {
  override def readFile(path: String): String = throw new NotImplementedError("readFile not available in OnlyPureFunctions.")
  override def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue] = throw new NotImplementedError("read_json not available in OnlyPureFunctions.")
  override def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError("write_json not available in OnlyPureFunctions.")
  override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = throw new NotImplementedError("size not available in OnlyPureFunctions.")
  override def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError("write_tsv not available in OnlyPureFunctions.")
  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError("stdout not available in OnlyPureFunctions.")
  override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedError("glob not available in OnlyPureFunctions.")
  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = throw new NotImplementedError("writeTempFile not available in OnlyPureFunctions.")
  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError("stderr not available in OnlyPureFunctions.")
}

trait PureFunctions { this: WdlStandardLibraryFunctions =>

  override def transpose(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    def extractExactlyOneArg: Try[WdlValue] = params.size match {
      case 1 => params.head
      case n => Failure(new IllegalArgumentException(s"Invalid number of parameters for engine function transpose: $n. Ensure transpose(x: Array[Array[X]]) takes exactly 1 parameters."))
    }

    case class ExpandedTwoDimensionalArray(innerType: WdlType, value: Seq[Seq[WdlValue]])
    def validateAndExpand(value: WdlValue): Try[ExpandedTwoDimensionalArray] = value match {
      case WdlArray(WdlArrayType(WdlArrayType(innerType)), array: Seq[WdlValue]) => expandWdlArray(array) map { ExpandedTwoDimensionalArray(innerType, _) }
      case array @ WdlArray(WdlArrayType(nonArrayType), _) => Failure(new IllegalArgumentException(s"Array must be two-dimensional to be transposed but given array of $nonArrayType"))
      case otherValue => Failure(new IllegalArgumentException(s"Function 'transpose' must be given a two-dimensional array but instead got ${otherValue.typeName}"))
    }

    def expandWdlArray(outerArray: Seq[WdlValue]): Try[Seq[Seq[WdlValue]]] = Try {
      outerArray map {
        case array: WdlArray => array.value
        case otherValue => throw new IllegalArgumentException(s"Function 'transpose' must be given a two-dimensional array but instead got WdlArray[${otherValue.typeName}]")
      }
    }

    def transpose(expandedTwoDimensionalArray: ExpandedTwoDimensionalArray): Try[WdlArray] = Try {
      val innerType = expandedTwoDimensionalArray.innerType
      val array = expandedTwoDimensionalArray.value
      WdlArray(WdlArrayType(WdlArrayType(innerType)), array.transpose map { WdlArray(WdlArrayType(innerType), _) })
    }

    extractExactlyOneArg.flatMap(validateAndExpand).flatMap(transpose)
  }

  override def range(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    def extractAndValidateArguments = params.size match {
      case 1 => validateArguments(params.head)
      case n => Failure(new IllegalArgumentException(s"Invalid number of parameters for engine function range: $n. Ensure range(x: WdlInteger) takes exactly 1 parameters."))
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
