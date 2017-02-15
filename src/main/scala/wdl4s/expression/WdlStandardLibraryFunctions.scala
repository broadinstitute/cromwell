package wdl4s.expression

import cats.instances.try_._
import cats.syntax.cartesian._
import wdl4s.types._
import wdl4s.values._
import wdl4s.{TsvSerializable, WdlExpressionException}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import WdlStandardLibraryFunctions.{crossProduct => stdLibCrossProduct, _}

trait WdlStandardLibraryFunctions extends WdlFunctions[WdlValue] {
  def fileContentsToString(path: String): String = readFile(path)
  def readFile(path: String): String
  def writeTempFile(path: String, prefix: String, suffix: String, content: String): String
  def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile]
  def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile]
  def glob(path: String, pattern: String): Seq[String]
  def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile]
  def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile]
  def size(params: Seq[Try[WdlValue]]): Try[WdlFloat]

  def read_objects(params: Seq[Try[WdlValue]]): Try[WdlArray] = extractObjects("read_objects", params) map { WdlArray(WdlArrayType(WdlObjectType), _) }
  def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = readContentsFromSingleFileParameter("read_string", params).map(s => WdlString(s.trim))
  def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue]
  def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] = read_string(params) map { s => WdlInteger(s.value.trim.toInt) }
  def read_float(params: Seq[Try[WdlValue]]): Try[WdlFloat] = read_string(params) map { s => WdlFloat(s.value.trim.toDouble) }

  def write_lines(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlArray])
  def write_map(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlMap])
  def write_object(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlObject])
  def write_objects(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv(params, classOf[WdlArray])

  def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      contents <- readContentsFromSingleFileParameter("read_lines", params)
      lines = contents.split("\n")
    } yield WdlArray(WdlArrayType(WdlStringType), lines map WdlString)
  }

  def read_map(params: Seq[Try[WdlValue]]): Try[WdlMap] = {
    for {
      contents <- readContentsFromSingleFileParameter("read_map", params)
      wdlMap <- WdlMap.fromTsv(contents, WdlMapType(WdlAnyType, WdlAnyType))
    } yield wdlMap
  }

  def read_object(params: Seq[Try[WdlValue]]): Try[WdlObject] = {
    extractObjects("read_object", params) map {
      case array if array.length == 1 => array.head
      case _ => throw new IllegalArgumentException("read_object yields an Object and thus can only read 2-rows TSV files. Try using read_objects instead.")
    }
  }

  def read_tsv(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      contents <- readContentsFromSingleFileParameter("read_tsv", params)
      wdlArray = WdlArray.fromTsv(contents)
    } yield wdlArray
  }

  def read_boolean(params: Seq[Try[WdlValue]]): Try[WdlBoolean] = {
    read_string(params) map { s => WdlBoolean(java.lang.Boolean.parseBoolean(s.value.trim.toLowerCase)) }
  }

  def glob(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    for {
      singleArgument <- extractSingleArgument("glob", params)
      globVal = singleArgument.valueString
      files = glob(globPath(globVal), globVal)
      wdlFiles = files map { WdlFile(_, isGlob = false) }
    } yield WdlArray(WdlArrayType(WdlFileType), wdlFiles)
  }

  def transpose(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
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

  def length(params: Seq[Try[WdlValue]]): Try[WdlInteger] = {
    def extractArguments: Try[WdlValue] = params.size match {
      case 1 => params.head
      case n => Failure(new IllegalArgumentException(s"Invalid number of parameters for engine function length(): $n. length() takes exactly 1 parameter."))
    }

    def arrayLength(value: WdlValue): Try[WdlInteger] = value match {
      case WdlArray(_, arrayValues) => Success(WdlInteger(arrayValues.length))
      case bad => Failure(new UnsupportedOperationException(s"length() expects one parameter of type Array but got one parameter of type ${bad.wdlType.toWdlString}"))
    }

    extractArguments flatMap arrayLength
  }

  def prefix(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    def extractTwoArguments: Try[(WdlValue, WdlValue)] = params.size match {
      case 2 => (params.head |@| params.tail.head) map { (_, _) }
      case n => Failure(new UnsupportedOperationException(s"prefix() expects two parameters but got $n"))
    }

    val makePrefixedString = (prefixString: WdlValue, elements: WdlValue) => (prefixString, elements) match {
      case (WdlString(p), WdlArray(WdlArrayType(etype), es)) if etype.isInstanceOf[WdlPrimitiveType] =>
        val result = es map { e => WdlString(p + e.valueString) }
        Success(WdlArray(WdlArrayType(WdlStringType), result))
      case (x, y) =>
        Failure(new UnsupportedOperationException(s"The function prefix expect arguments (String, Array[X]) where X is a primitive type, but got (${prefixString.wdlType.toWdlString}, ${elements.wdlType.toWdlString})"))
    }

    extractTwoArguments flatMap makePrefixedString.tupled
  }

  def range(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
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

  def sub(params: Seq[Try[WdlValue]]): Try[WdlString] = {
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

  private val SelectFirstEmptyInput = Failure(new IllegalArgumentException("select_first failed. The input array was empty."))
  def select_first(params: Seq[Try[WdlValue]]): Try[WdlValue] = extractSingleArgument("select_first", params) flatMap {
    case WdlArray(WdlArrayType(WdlOptionalType(memberType)), arrayValue) =>
      if (arrayValue.isEmpty) SelectFirstEmptyInput else (arrayValue collectFirst {
        case WdlOptionalValue(_, Some(wdlValue)) => wdlValue
        case wdlValue if memberType.isCoerceableFrom(wdlValue.wdlType) =>
          memberType.coerceRawValue(wdlValue).get}).map(Success(_)).getOrElse(Failure(new IllegalArgumentException("select_first failed. All provided values were empty.")))
    case allValid @ WdlArray(WdlArrayType(_), arrayValue) => if (arrayValue.isEmpty) SelectFirstEmptyInput else Success(arrayValue.head)
    case other => Failure(new IllegalArgumentException(s"select_first must take an array but got ${other.wdlType.toWdlString}: ${other.toWdlString}"))
  }

  def select_all(params: Seq[Try[WdlValue]]): Try[WdlArray] = extractSingleArgument("select_all", params) flatMap {
    case WdlArray(WdlArrayType(WdlOptionalType(memberType)), arrayValue) =>
      Success(WdlArray(WdlArrayType(memberType), arrayValue collect {
        case WdlOptionalValue(_, Some(wdlValue)) => wdlValue
        case wdlValue if memberType.isCoerceableFrom(wdlValue.wdlType) => memberType.coerceRawValue(wdlValue).get
      }))
    case allValid @ WdlArray(WdlArrayType(_), _) => Success(allValid)
    case other => Failure(new IllegalArgumentException("select_all must take an array but got: " + other.toWdlString))
  }

  def defined(params: Seq[Try[WdlValue]]): Try[WdlBoolean] = extractSingleArgument("defined", params) map {
    case WdlOptionalValue(_, None) => WdlBoolean(false)
    case _ => WdlBoolean(true)
  }

  def zip(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    val badArgsFailure = Failure(new IllegalArgumentException(s"Invalid parameters for engine function zip: $params. Requires exactly two evaluated array values of equal length."))

    for {
      values <- extractTwoParams(params, badArgsFailure)
      (left, right) <- assertEquallySizedArrays(values, badArgsFailure)
      leftType = left.wdlType.memberType
      rightType = right.wdlType.memberType
      zipped = left.value.zip(right.value) map { case (l,r) => WdlPair(l, r) }
    } yield WdlArray(WdlArrayType(WdlPairType(leftType, rightType)), zipped)
  }

  def cross(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    val badArgsFailure = Failure(new IllegalArgumentException(s"Invalid parameters for engine function cross: $params. Requires exactly two evaluated array values of equal length."))

    for {
      values <- extractTwoParams(params, badArgsFailure)
      (left, right) <- assertArrays(values, badArgsFailure)
      leftType = left.wdlType.memberType
      rightType = right.wdlType.memberType
      crossed = stdLibCrossProduct(left.value, right.value) map { case (l,r) => WdlPair(l, r) }
    } yield WdlArray(WdlArrayType(WdlPairType(leftType, rightType)), crossed)
  }

  /**
    * Asserts that the parameter list contains a single parameter which will be interpreted
    * as a File and attempts to read the contents of that file and returns back the contents
    * as a String
    */
  private def readContentsFromSingleFileParameter(functionName: String, params: Seq[Try[WdlValue]]): Try[String] = {
    for {
      singleArgument <- extractSingleArgument(functionName, params)
      string = fileContentsToString(singleArgument.valueString)
    } yield string
  }

  private def extractObjects(functionName: String, params: Seq[Try[WdlValue]]): Try[Array[WdlObject]] = for {
    contents <- readContentsFromSingleFileParameter(functionName, params)
    wdlObjects <- WdlObject.fromTsv(contents)
  } yield wdlObjects

  private def writeContent(baseName: String, content: String): Try[WdlFile] = {
    Try(WdlFile(writeTempFile(tempFilePath, s"$baseName.", ".tmp", content)))
  }

  private def writeToTsv(params: Seq[Try[WdlValue]], wdlClass: Class[_ <: WdlValue with TsvSerializable]) = {
    for {
      singleArgument <- extractSingleArgument("writeToTsv", params)
      downcast <- Try(wdlClass.cast(singleArgument))
      tsvSerialized <- downcast.tsvSerialize
      file <- writeContent(wdlClass.getSimpleName.toLowerCase, tsvSerialized)
    } yield file
  }
}

object WdlStandardLibraryFunctions {
  def crossProduct[A, B](as: Seq[A], bs: Seq[B]): Seq[(A, B)] = for {
    a <- as
    b <- bs
  } yield (a, b)

  def extractTwoParams[A](params: Seq[Try[A]], badArgsFailure: Failure[Nothing]): Try[(A, A)] = {
    if (params.size != 2) { badArgsFailure }
    else for {
      left <- params.head
      right <- params(1)
    } yield (left, right)
  }

  def assertEquallySizedArrays[A](values: (WdlValue, WdlValue), badArgsFailure: Failure[Nothing] ): Try[(WdlArray, WdlArray)] = values match {
    case (leftArray: WdlArray, rightArray: WdlArray) if leftArray.value.size == rightArray.value.size => Success((leftArray, rightArray))
    case _ => badArgsFailure
  }

  def assertArrays(values: (WdlValue, WdlValue), badArgsFailure: Failure[Nothing] ): Try[(WdlArray, WdlArray)] = values match {
    case (leftArray: WdlArray, rightArray: WdlArray) => Success((leftArray, rightArray))
    case _ => badArgsFailure
  }
}

trait PureStandardLibraryFunctionsLike extends WdlStandardLibraryFunctions {

  def className = this.getClass.getCanonicalName

  override def readFile(path: String): String = throw new NotImplementedError(s"readFile not available in $className.")
  override def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue] = throw new NotImplementedError(s"read_json not available in $className.")
  override def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError(s"write_json not available in $className.")
  override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = throw new NotImplementedError(s"size not available in $className.")
  override def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError(s"write_tsv not available in $className.")
  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError(s"stdout not available in $className.")
  override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedError(s"glob not available in $className.")
  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = throw new NotImplementedError(s"writeTempFile not available in $className.")
  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = throw new NotImplementedError(s"stderr not available in $className.")
}

case object PureStandardLibraryFunctions extends PureStandardLibraryFunctionsLike

class WdlStandardLibraryFunctionsType extends WdlFunctions[WdlType] {
  def stdout(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def stderr(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def read_lines(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlStringType))
  def read_tsv(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlArrayType(WdlStringType)))
  def read_map(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlMapType(WdlStringType, WdlStringType))
  def read_object(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlObjectType)
  def read_objects(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlObjectType))
  def read_json(params: Seq[Try[WdlType]]): Try[WdlType] = Failure(new WdlExpressionException("Return type of read_json() can't be known statically"))
  def read_int(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlIntegerType)
  def read_string(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlStringType)
  def read_float(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFloatType)
  def read_boolean(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlBooleanType)
  def write_lines(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def write_tsv(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def write_map(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def write_object(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def write_objects(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def write_json(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFileType)
  def glob(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlFileType))
  def size(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlFloatType)
  def length(params: Seq[Try[WdlType]]): Try[WdlType] = params.toList match {
    case Success(WdlArrayType(_)) :: Nil => Success(WdlIntegerType)
    case _ =>
      val badArgs = params.mkString(", ")
      Failure(new Exception(s"Unexpected arguments to function `length`.  `length` takes a parameter of type Array but got: $badArgs"))
  }
  def prefix(params: Seq[Try[WdlType]]): Try[WdlType] = params.toList match {
    case Success(WdlStringType) :: Success(WdlArrayType(_: WdlPrimitiveType)) :: Nil => Success(WdlArrayType(WdlStringType))
    case _ =>
      val badArgs = params.mkString(", ")
      Failure(new Exception(s"Unexpected arguments to function `prefix`.  `prefix` takes parameters of type String and Array[<primitive>] but got: $badArgs"))
  }
  def sub(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlStringType)
  def range(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlArrayType(WdlIntegerType))
  def transpose(params: Seq[Try[WdlType]]): Try[WdlType] = params.toList match {
    case Success(t @ WdlArrayType(WdlArrayType(wdlType))) :: Nil => Success(t)
    case _ => Failure(new Exception(s"Unexpected transpose target: $params"))
  }
  def select_first(params: Seq[Try[WdlType]]): Try[WdlType] = extractSingleArgument("select_first", params) flatMap {
    case WdlArrayType(WdlOptionalType(innerType)) => Success(innerType)
    case other => Failure(new IllegalArgumentException(s"select_first failed. It expects an array of optional values but got ${other.toWdlString}."))
  }
  def select_all(params: Seq[Try[WdlType]]): Try[WdlType] = extractSingleArgument("select_all", params) flatMap {
    case WdlArrayType(WdlOptionalType(innerType)) => Success(WdlArrayType(innerType))
    case other => Failure(new IllegalArgumentException(s"select_all failed. It expects an array of optional values but got ${other.toWdlString}."))
  }
  def defined(params: Seq[Try[WdlType]]): Try[WdlType] = extractSingleArgument("defined", params).map(_ => WdlBooleanType)
  def zip(params: Seq[Try[WdlType]]): Try[WdlType] = {
    val badArgsFailure = Failure(new Exception(s"Unexpected zip parameters: $params"))
    WdlStandardLibraryFunctions.extractTwoParams(params, badArgsFailure) flatMap {
      case (arrayType1: WdlArrayType, arrayType2: WdlArrayType) => Success(WdlArrayType(WdlPairType(arrayType1.memberType, arrayType2.memberType)))
      case _ => badArgsFailure
    }
  }
  def cross(params: Seq[Try[WdlType]]): Try[WdlType] = zip(params)
}

case object NoFunctions extends WdlStandardLibraryFunctions {
  override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedError()
  override def readFile(path: String): String = throw new NotImplementedError()
  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = throw new NotImplementedError()
  override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError())
  override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError())
  override def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue] = Failure(new NotImplementedError())
  override def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError())
  override def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError())
  override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = Failure(new NotImplementedError())
  override def length(params: Seq[Try[WdlValue]]): Try[WdlInteger] = Failure(new NotImplementedError())
  override def sub(params: Seq[Try[WdlValue]]): Try[WdlString] = Failure(new NotImplementedError())
  override def range(params: Seq[Try[WdlValue]]): Try[WdlArray] = Failure(new NotImplementedError())
  override def transpose(params: Seq[Try[WdlValue]]): Try[WdlArray] = Failure(new NotImplementedError())
  override def select_first(params: Seq[Try[WdlValue]]): Try[WdlArray] = Failure(new NotImplementedError())
  override def select_all(params: Seq[Try[WdlValue]]): Try[WdlArray] = Failure(new NotImplementedError())
  override def zip(params: Seq[Try[WdlValue]]): Try[WdlArray] = Failure(new NotImplementedError())
  override def cross(params: Seq[Try[WdlValue]]): Try[WdlArray] = Failure(new NotImplementedError())
}
