package wdl.expression

import cats.instances.try_._
import cats.syntax.apply._
import common.exception.AggregatedException
import common.util.TryUtil
import wdl.expression.WdlStandardLibraryFunctions.{crossProduct => stdLibCrossProduct, _}
import wom.TsvSerializable
import wom.expression.IoFunctionSet
import wom.types._
import wom.values.WomArray.WomArrayLike
import wom.values._
import spray.json._

import scala.concurrent.Await
import scala.concurrent.duration.Duration
import scala.util.{Failure, Success, Try}

trait WdlStandardLibraryFunctions extends WdlFunctions[WomValue] {
  def readFile(path: String): String

  def writeFile(path: String, content: String): Try[WomFile]

  def stdout(params: Seq[Try[WomValue]]): Try[WomFile]

  def stderr(params: Seq[Try[WomValue]]): Try[WomFile]

  def size(params: Seq[Try[WomValue]]): Try[WomFloat]

  private def writeContent(baseName: String, content: String): Try[WomFile] = writeFile(s"${baseName}_${content.md5Sum}.tmp", content)

  private def writeToTsv[A <: WomValue with TsvSerializable](functionName: String, params: Seq[Try[WomValue]], defaultIfOptionalEmpty: A): Try[WomFile] = {
    val wdlClass = defaultIfOptionalEmpty.getClass
    def castOrDefault(womValue: WomValue): A = womValue match {
      case WomOptionalValue(_, None) => defaultIfOptionalEmpty
      case WomOptionalValue(_, Some(v)) => wdlClass.cast(v)
      case _ => wdlClass.cast(womValue)
    }

    for {
      singleArgument <- extractSingleArgument(functionName, params)
      downcast <- Try(castOrDefault(singleArgument))
      tsvSerialized <- downcast.tsvSerialize
      file <- writeContent(functionName, tsvSerialized)
    } yield file
  }

  def read_objects(params: Seq[Try[WomValue]]): Try[WomArray] = extractObjects("read_objects", params) map { WomArray(WomArrayType(WomObjectType), _) }
  def read_string(params: Seq[Try[WomValue]]): Try[WomString] = readContentsFromSingleFileParameter("read_string", params).map(s => WomString(s.trim))
  def read_json(params: Seq[Try[WomValue]]): Try[WomValue] = readContentsFromSingleFileParameter("read_json", params).map(_.parseJson).flatMap(WomAnyType.coerceRawValue)
  def read_int(params: Seq[Try[WomValue]]): Try[WomInteger] = read_string(params) map { s => WomInteger(s.value.trim.toInt) }
  def read_float(params: Seq[Try[WomValue]]): Try[WomFloat] = read_string(params) map { s => WomFloat(s.value.trim.toDouble) }

  def write_lines(params: Seq[Try[WomValue]]): Try[WomFile] = writeToTsv("write_lines", params, WomArray(WomArrayType(WomStringType), List.empty[WomValue]))
  def write_map(params: Seq[Try[WomValue]]): Try[WomFile] = writeToTsv("write_map", params, WomMap(WomMapType(WomStringType, WomStringType), Map.empty[WomValue, WomValue]))
  def write_object(params: Seq[Try[WomValue]]): Try[WomFile] = writeToTsv("write_object", params, WomObject(Map.empty[String, WomValue]))
  def write_objects(params: Seq[Try[WomValue]]): Try[WomFile] = writeToTsv("write_objects", params, WomArray(WomArrayType(WomObjectType), List.empty[WomObject]))
  def write_tsv(params: Seq[Try[WomValue]]): Try[WomFile] = writeToTsv("write_tsv", params, WomArray(WomArrayType(WomStringType), List.empty[WomValue]))
  def write_json(params: Seq[Try[WomValue]]): Try[WomFile] = for {
    value <- extractSingleArgument("write_json", params)
    jsonContent = valueToJson(value)
    written <- writeContent("write_json", jsonContent.compactPrint)
  } yield written

  private def valueToJson(womValue: WomValue): JsValue = womValue match {
    case WomInteger(i) => JsNumber(i)
    case WomFloat(f) => JsNumber(f)
    case WomString(s) => JsString(s)
    case WomBoolean(b) => JsBoolean(b)
    case f: WomFile => JsString(f.value)
    case WomPair(left, right) => JsObject(Map("left" -> valueToJson(left), "right" -> valueToJson(right)))
    case WomArray(_, values) => JsArray(values.map(valueToJson).toVector)
    case WomMap(_, value) => JsObject(value map { case (k, v) => k.valueString -> valueToJson(v) })
    case o: WomObjectLike => JsObject(o.values map { case (k, v) => k -> valueToJson(v) })
    case opt: WomOptionalValue => opt.value match {
      case Some(inner) => valueToJson(inner)
      case None => JsNull
    }
  }

  def read_lines(params: Seq[Try[WomValue]]): Try[WomArray] = {
    for {
      contents <- readContentsFromSingleFileParameter("read_lines", params)
      lines = contents.split("\n")
    } yield WomArray(WomArrayType(WomStringType), lines map WomString)
  }

  def read_map(params: Seq[Try[WomValue]]): Try[WomMap] = {
    for {
      contents <- readContentsFromSingleFileParameter("read_map", params)
      wdlMap <- WomMap.fromTsv(contents, WomMapType(WomAnyType, WomAnyType))
    } yield wdlMap
  }

  def read_object(params: Seq[Try[WomValue]]): Try[WomObject] = {
    extractObjects("read_object", params) map {
      case array if array.length == 1 => array.head
      case _ => throw new IllegalArgumentException("read_object yields an Object and thus can only read 2-rows TSV files. Try using read_objects instead.")
    }
  }

  def read_tsv(params: Seq[Try[WomValue]]): Try[WomArray] = {
    for {
      contents <- readContentsFromSingleFileParameter("read_tsv", params)
      wdlArray = WomArray.fromTsv(contents)
    } yield wdlArray
  }

  def read_boolean(params: Seq[Try[WomValue]]): Try[WomBoolean] = {
    read_string(params) map { s => WomBoolean(java.lang.Boolean.parseBoolean(s.value.trim.toLowerCase)) }
  }

  def globHelper(pattern: String): Seq[String]

  final def glob(params: Seq[Try[WomValue]]): Try[WomArray] = {
    for {
      pattern <- extractSingleArgument("glob", params)
      womString <- WomStringType.coerceRawValue(pattern)
      patternString = womString.valueString
      filePaths <- Try(globHelper(patternString))
    } yield WomArray(WomArrayType(WomSingleFileType), filePaths.map(WomSingleFile(_)))
  }

  def basename(params: Seq[Try[WomValue]]): Try[WomString] = {

    val arguments: Try[(WomValue, WomValue)] = params.toList match {
      case Success(f) :: Nil => Success((f, WomString("")))
      case Success(f) :: Success(s) :: Nil => Success((f, s))
      case s if s.lengthCompare(2) > 0 || s.lengthCompare(1) < 0 =>
        val message = s"Bad number of arguments to basename(filename, suffixToStrip = ''): ${params.size}"
        Failure(new IllegalArgumentException(message))
      case _ =>
        val failures = params collect {
          case Failure(e) => e
        }
        Failure(AggregatedException("Failures evaluating basename parameters", failures))
    }

    for {
      extractedArgs <- arguments
      fileNameAsString <- WomStringType.coerceRawValue(extractedArgs._1)
      suffixAsString <- WomStringType.coerceRawValue(extractedArgs._2)
      basename = fileNameAsString.valueString.split('/').last
      suffixless = basename.stripSuffix(suffixAsString.valueString)
    } yield WomString(suffixless)
  }

  def floor(params: Seq[Try[WomValue]]): Try[WomInteger] = {
    extractSingleArgument("floor", params) flatMap { f => WomFloatType.coerceRawValue(f) } map { f => WomInteger(Math.floor(f.asInstanceOf[WomFloat].value).toInt) }
  }

  def round(params: Seq[Try[WomValue]]): Try[WomInteger] = {
    extractSingleArgument("round", params) flatMap { f => WomFloatType.coerceRawValue(f) } map { f => WomInteger(Math.round(f.asInstanceOf[WomFloat].value).toInt) }
  }

  def ceil(params: Seq[Try[WomValue]]): Try[WomInteger] = {
    extractSingleArgument("ceil", params) flatMap { f => WomFloatType.coerceRawValue(f) } map { f => WomInteger(Math.ceil(f.asInstanceOf[WomFloat].value).toInt) }
  }

  def transpose(params: Seq[Try[WomValue]]): Try[WomArray] = {
    def extractExactlyOneArg: Try[WomValue] = params.size match {
      case 1 => params.head
      case n => Failure(new IllegalArgumentException(s"Invalid number of parameters for engine function transpose: $n. Ensure transpose(x: Array[Array[X]]) takes exactly 1 parameters."))
    }

    case class ExpandedTwoDimensionalArray(innerType: WomType, value: Seq[Seq[WomValue]])
    def validateAndExpand(value: WomValue): Try[ExpandedTwoDimensionalArray] = value match {
      case WomArray(WomArrayType(WomArrayType(innerType)), array: Seq[WomValue]) => expandWdlArray(array) map { ExpandedTwoDimensionalArray(innerType, _) }
      case WomArray(WomArrayType(nonArrayType), _) => Failure(new IllegalArgumentException(s"Array must be two-dimensional to be transposed but given array of $nonArrayType"))
      case otherValue => Failure(new IllegalArgumentException(s"Function 'transpose' must be given a two-dimensional array but instead got ${otherValue.typeName}"))
    }

    def expandWdlArray(outerArray: Seq[WomValue]): Try[Seq[Seq[WomValue]]] = Try {
      outerArray map {
        case array: WomArray => array.value
        case otherValue => throw new IllegalArgumentException(s"Function 'transpose' must be given a two-dimensional array but instead got WdlArray[${otherValue.typeName}]")
      }
    }

    def transpose(expandedTwoDimensionalArray: ExpandedTwoDimensionalArray): Try[WomArray] = Try {
      val innerType = expandedTwoDimensionalArray.innerType
      val array = expandedTwoDimensionalArray.value
      WomArray(WomArrayType(WomArrayType(innerType)), array.transpose map { WomArray(WomArrayType(innerType), _) })
    }

    extractExactlyOneArg.flatMap(validateAndExpand).flatMap(transpose)
  }

  def length(params: Seq[Try[WomValue]]): Try[WomInteger] = {
    def extractArguments: Try[WomValue] = params.size match {
      case 1 => params.head
      case n => Failure(new IllegalArgumentException(s"Invalid number of parameters for engine function length(): $n. length() takes exactly 1 parameter."))
    }

    def arrayLength(value: WomValue): Try[WomInteger] = value match {
      case WomArray(_, arrayValues) => Success(WomInteger(arrayValues.length))
      case bad => Failure(new UnsupportedOperationException(s"length() expects one parameter of type Array but got one parameter of type ${bad.womType.toDisplayString}"))
    }

    extractArguments flatMap arrayLength
  }

  def flatten(params: Seq[Try[WomValue]]): Try[WomValue] = {
    def getFlatValues(v: WomValue): Try[Seq[WomValue]] = v match {
      case WomArrayLike(WomArray(_, values)) => Success(values.toList)
      case other => Failure(new IllegalArgumentException(s"Invalid argument to flatten: ${other.womType.toDisplayString}, flatten requires an Array[Array[_]]"))
    }

    val arg: Try[WomValue] = extractSingleArgument("flatten", params)
    arg flatMap {
      case WomArray(WomArrayType(WomArrayType(elemType)), arrayValues) =>
        val llt: Try[Seq[Seq[WomValue]]] = TryUtil.sequence(arrayValues.map(getFlatValues))
        llt.map(ll => WomArray(WomArrayType(elemType), ll.flatten))
      case bad =>
        Failure(new UnsupportedOperationException(s"flatten() expects one parameter of type Array[Array[T]] but got one parameter of type ${bad.womType.toDisplayString}"))
    }
  }

  def prefix(params: Seq[Try[WomValue]]): Try[WomArray] = {
    def extractTwoArguments: Try[(WomValue, WomValue)] = params.size match {
      case 2 => (params.head, params.tail.head).tupled
      case n => Failure(new UnsupportedOperationException(s"prefix() expects two parameters but got $n"))
    }

    val makePrefixedString = (prefixString: WomValue, elements: WomValue) => (prefixString, elements) match {
      case (WomString(p), WomArray(WomArrayType(etype), es)) if etype.isInstanceOf[WomPrimitiveType] =>
        val result = es map { e => WomString(p + e.valueString) }
        Success(WomArray(WomArrayType(WomStringType), result))
      case (_, _) =>
        Failure(new UnsupportedOperationException(s"The function prefix expect arguments (String, Array[X]) where X is a primitive type, but got (${prefixString.womType.toDisplayString}, ${elements.womType.toDisplayString})"))
    }

    extractTwoArguments flatMap makePrefixedString.tupled
  }

  def range(params: Seq[Try[WomValue]]): Try[WomArray] = {
    def extractAndValidateArguments = params.size match {
      case 1 => validateArguments(params.head)
      case n => Failure(new IllegalArgumentException(s"Invalid number of parameters for engine function range: $n. Ensure range(x: WdlInteger) takes exactly 1 parameters."))
    }

    def validateArguments(value: Try[WomValue]) = value match {
      case Success(intValue: WomValue) if WomIntegerType.isCoerceableFrom(intValue.womType) =>
        Integer.valueOf(intValue.valueString) match {
          case i if i >= 0 => Success(i)
          case n => Failure(new IllegalArgumentException(s"Parameter to seq must be greater than or equal to 0 (but got $n)"))
        }
      case _ => Failure(new IllegalArgumentException(s"Invalid parameter for engine function seq: $value."))
    }

    extractAndValidateArguments map { intValue => WomArray(WomArrayType(WomIntegerType), (0 until intValue).map(WomInteger)) }
  }

  def sub(params: Seq[Try[WomValue]]): Try[WomString] = {
    def extractArguments = params.size match {
      case 3 => Success((params.head, params(1), params(2)))
      case n => Failure(new IllegalArgumentException(s"Invalid number of parameters for engine function sub: $n. sub takes exactly 3 parameters."))
    }

    def validateArguments(values: (Try[WomValue], Try[WomValue], Try[WomValue])) = values match {
      case (Success(strValue), Success(WomString(pattern)), Success(replaceValue))
        if WomStringType.isCoerceableFrom(strValue.womType) &&
          WomStringType.isCoerceableFrom(replaceValue.womType) =>
        Success((strValue.valueString, pattern, replaceValue.valueString))
      case _ => Failure(new IllegalArgumentException(s"Invalid parameters for engine function sub: $values."))
    }

    for {
      args <- extractArguments
      (str, pattern, replace) <- validateArguments(args)
    } yield WomString(pattern.r.replaceAllIn(str, replace))
  }

  private val SelectFirstEmptyInput = Failure(new IllegalArgumentException("select_first failed. The input array was empty."))
  def select_first(params: Seq[Try[WomValue]]): Try[WomValue] = extractSingleArgument("select_first", params) flatMap {
    case WomArray(WomArrayType(WomOptionalType(memberType)), arrayValue) =>
      if (arrayValue.isEmpty) SelectFirstEmptyInput else (arrayValue collectFirst {
        case WomOptionalValue(_, Some(womValue)) => womValue
        case womValue if memberType.isCoerceableFrom(womValue.womType) =>
          memberType.coerceRawValue(womValue).get}).map(Success(_)).getOrElse(Failure(new IllegalArgumentException("select_first failed. All provided values were empty.")))
    case WomArray(WomArrayType(_), arrayValue) => if (arrayValue.isEmpty) SelectFirstEmptyInput else Success(arrayValue.head)
    case other => Failure(new IllegalArgumentException(s"select_first must take an array but got ${other.womType.toDisplayString}: ${other.toWomString}"))
  }

  def select_all(params: Seq[Try[WomValue]]): Try[WomArray] = extractSingleArgument("select_all", params) flatMap {
    case WomArray(WomArrayType(WomOptionalType(memberType)), arrayValue) =>
      Success(WomArray(WomArrayType(memberType), arrayValue collect {
        case WomOptionalValue(_, Some(womValue)) => womValue
        case womValue if memberType.isCoerceableFrom(womValue.womType) => memberType.coerceRawValue(womValue).get
      }))
    case allValid @ WomArray(WomArrayType(_), _) => Success(allValid)
    case other => Failure(new IllegalArgumentException("select_all must take an array but got: " + other.toWomString))
  }

  def defined(params: Seq[Try[WomValue]]): Try[WomBoolean] = extractSingleArgument("defined", params) map {
    case WomOptionalValue(_, None) => WomBoolean(false)
    case _ => WomBoolean(true)
  }

  def zip(params: Seq[Try[WomValue]]): Try[WomArray] = {
    val badArgsFailure = Failure(new IllegalArgumentException(s"Invalid parameters for engine function zip: $params. Requires exactly two evaluated array values of equal length."))

    for {
      values <- extractTwoParams(params, badArgsFailure)
      (left, right) <- assertEquallySizedArrays(values, badArgsFailure)
      leftType = left.womType.memberType
      rightType = right.womType.memberType
      zipped = left.value.zip(right.value) map { case (l,r) => WomPair(l, r) }
    } yield WomArray(WomArrayType(WomPairType(leftType, rightType)), zipped)
  }

  def cross(params: Seq[Try[WomValue]]): Try[WomArray] = {
    val badArgsFailure = Failure(new IllegalArgumentException(s"Invalid parameters for engine function cross: $params. Requires exactly two evaluated array values of equal length."))

    for {
      values <- extractTwoParams(params, badArgsFailure)
      (left, right) <- assertArrays(values, badArgsFailure)
      leftType = left.womType.memberType
      rightType = right.womType.memberType
      crossed = stdLibCrossProduct(left.value, right.value) map { case (l,r) => WomPair(l, r) }
    } yield WomArray(WomArrayType(WomPairType(leftType, rightType)), crossed)
  }

  /**
    * Asserts that the parameter list contains a single parameter which will be interpreted
    * as a File and attempts to read the contents of that file and returns back the contents
    * as a String
    */
  private def readContentsFromSingleFileParameter(functionName: String, params: Seq[Try[WomValue]]): Try[String] = {
    for {
      singleArgument <- extractSingleArgument(functionName, params)
      string = readFile(singleArgument.valueString)
    } yield string
  }

  private def extractObjects(functionName: String, params: Seq[Try[WomValue]]): Try[Array[WomObject]] = for {
    contents <- readContentsFromSingleFileParameter(functionName, params)
    wdlObjects <- WomObject.fromTsv(contents)
  } yield wdlObjects
}

object WdlStandardLibraryFunctions {
  def fromIoFunctionSet(ioFunctionSet: IoFunctionSet) = new WdlStandardLibraryFunctions {
    override def readFile(path: String): String = Await.result(ioFunctionSet.readFile(path, None, failOnOverflow = false), Duration.Inf)

    override def writeFile(path: String, content: String): Try[WomFile] = Try(Await.result(ioFunctionSet.writeFile(path, content), Duration.Inf))

    override def stdout(params: Seq[Try[WomValue]]): Try[WomFile] = ioFunctionSet.stdout(params)

    override def stderr(params: Seq[Try[WomValue]]): Try[WomFile] = ioFunctionSet.stderr(params)

    override def globHelper(pattern: String): Seq[String] = Await.result(ioFunctionSet.glob(pattern), Duration.Inf)

    override def size(params: Seq[Try[WomValue]]): Try[WomFloat] = Try(Await.result(ioFunctionSet.size(params), Duration.Inf))
  }

  def crossProduct[A, B](as: Seq[A], bs: Seq[B]): Seq[(A, B)] = for {
    a <- as
    b <- bs
  } yield (a, b)

  def extractTwoParams[A](params: Seq[Try[A]], badArgsFailure: Failure[Nothing]): Try[(A, A)] = {
    if (params.lengthCompare(2) != 0) { badArgsFailure }
    else for {
      left <- params.head
      right <- params(1)
    } yield (left, right)
  }

  def assertEquallySizedArrays[A](values: (WomValue, WomValue), badArgsFailure: Failure[Nothing] ): Try[(WomArray, WomArray)] = values match {
    case (leftArray: WomArray, rightArray: WomArray) if leftArray.value.lengthCompare(rightArray.value.size) == 0 => Success((leftArray, rightArray))
    case _ => badArgsFailure
  }

  def assertArrays(values: (WomValue, WomValue), badArgsFailure: Failure[Nothing] ): Try[(WomArray, WomArray)] = values match {
    case (leftArray: WomArray, rightArray: WomArray) => Success((leftArray, rightArray))
    case _ => badArgsFailure
  }
}

trait PureStandardLibraryFunctionsLike extends WdlStandardLibraryFunctions {

  def className = this.getClass.getCanonicalName

  override def readFile(path: String): String = throw new NotImplementedError(s"readFile not available in $className.")
  override def writeFile(path: String, content: String): Try[WomFile] = throw new NotImplementedError(s"writeFile not available in $className.")
  override def read_json(params: Seq[Try[WomValue]]): Try[WomValue] = Failure(new NotImplementedError(s"read_json not available in $className."))
  override def write_json(params: Seq[Try[WomValue]]): Try[WomFile] = Failure(new NotImplementedError(s"write_json not available in $className."))
  override def size(params: Seq[Try[WomValue]]): Try[WomFloat] = Failure(new NotImplementedError(s"size not available in $className."))
  override def write_tsv(params: Seq[Try[WomValue]]): Try[WomFile] = Failure(new NotImplementedError(s"write_tsv not available in $className."))
  override def stdout(params: Seq[Try[WomValue]]): Try[WomFile] = Failure(new NotImplementedError(s"stdout not available in $className."))
  override def globHelper(pattern: String): Seq[String] = throw new NotImplementedError(s"glob not available in $className.")
  override def stderr(params: Seq[Try[WomValue]]): Try[WomFile] = Failure(new NotImplementedError(s"stderr not available in $className."))
}

case object PureStandardLibraryFunctions extends PureStandardLibraryFunctionsLike

class WdlStandardLibraryFunctionsType extends WdlFunctions[WomType] {
  def stdout(params: Seq[Try[WomType]]): Try[WomType] = Success(WomSingleFileType)
  def stderr(params: Seq[Try[WomType]]): Try[WomType] = Success(WomSingleFileType)
  def read_lines(params: Seq[Try[WomType]]): Try[WomType] = Success(WomArrayType(WomStringType))
  def read_tsv(params: Seq[Try[WomType]]): Try[WomType] = Success(WomArrayType(WomArrayType(WomStringType)))
  def read_map(params: Seq[Try[WomType]]): Try[WomType] = Success(WomMapType(WomStringType, WomStringType))
  def read_object(params: Seq[Try[WomType]]): Try[WomType] = Success(WomObjectType)
  def read_objects(params: Seq[Try[WomType]]): Try[WomType] = Success(WomArrayType(WomObjectType))
  def read_json(params: Seq[Try[WomType]]): Try[WomType] = Success(WomObjectType)
  def read_int(params: Seq[Try[WomType]]): Try[WomType] = Success(WomIntegerType)
  def read_string(params: Seq[Try[WomType]]): Try[WomType] = Success(WomStringType)
  def read_float(params: Seq[Try[WomType]]): Try[WomType] = Success(WomFloatType)
  def read_boolean(params: Seq[Try[WomType]]): Try[WomType] = Success(WomBooleanType)
  def write_lines(params: Seq[Try[WomType]]): Try[WomType] = Success(WomSingleFileType)
  def write_tsv(params: Seq[Try[WomType]]): Try[WomType] = Success(WomSingleFileType)
  def write_map(params: Seq[Try[WomType]]): Try[WomType] = Success(WomSingleFileType)
  def write_object(params: Seq[Try[WomType]]): Try[WomType] = Success(WomSingleFileType)
  def write_objects(params: Seq[Try[WomType]]): Try[WomType] = Success(WomSingleFileType)
  def write_json(params: Seq[Try[WomType]]): Try[WomType] = Success(WomSingleFileType)
  def glob(params: Seq[Try[WomType]]): Try[WomType] = Success(WomArrayType(WomSingleFileType))
  def size(params: Seq[Try[WomType]]): Try[WomType] = {
    def isGoodFirstSizeParam(womType: WomType): Boolean = womType match {
      case f if WomSingleFileType.isCoerceableFrom(f) => true
      case WomOptionalType(o) => isGoodFirstSizeParam(o)
      case _ => false
    }

    params.toList match {
      case Success(f) :: Nil if isGoodFirstSizeParam(f) => Success(WomFloatType)
      case Success(f) :: Success(WomStringType) :: Nil if isGoodFirstSizeParam(f) => Success(WomFloatType)
      case other => Failure(new Exception(s"Unexpected arguments to function `size`. Expected 'size(file: File [, unit: String])' but got 'size(${other.map(_.map(_.toDisplayString)).mkString(", ")})'"))
    }
  }
  def length(params: Seq[Try[WomType]]): Try[WomType] = params.toList match {
    case Success(WomArrayType(_)) :: Nil => Success(WomIntegerType)
    case _ =>
      val badArgs = params.mkString(", ")
      Failure(new Exception(s"Unexpected arguments to function `length`. `length` takes a parameter of type Array but got: $badArgs"))
  }
  def prefix(params: Seq[Try[WomType]]): Try[WomType] = params.toList match {
    case Success(WomStringType) :: Success(WomArrayType(_: WomPrimitiveType)) :: Nil => Success(WomArrayType(WomStringType))
    case _ =>
      val badArgs = params.mkString(", ")
      Failure(new Exception(s"Unexpected arguments to function `prefix`.  `prefix` takes parameters of type String and Array[<primitive>] but got: $badArgs"))
  }
  def sub(params: Seq[Try[WomType]]): Try[WomType] = Success(WomStringType)
  def range(params: Seq[Try[WomType]]): Try[WomType] = Success(WomArrayType(WomIntegerType))
  def floor(params: Seq[Try[WomType]]): Try[WomType] = params.toList match {
    case Success(ftype) :: Nil if WomFloatType.isCoerceableFrom(ftype) => Success(WomIntegerType)
    case _ => Failure(new Exception(s"Unexpected 'floor' arguments: $params (expects a single Float argument)"))
  }
  def round(params: Seq[Try[WomType]]): Try[WomType] = params.toList match {
    case Success(ftype) :: Nil if WomFloatType.isCoerceableFrom(ftype) => Success(WomIntegerType)
    case _ => Failure(new Exception(s"Unexpected 'round' arguments: $params (expects a single Float argument)"))
  }
  def ceil(params: Seq[Try[WomType]]): Try[WomType] = params.toList match {
    case Success(ftype) :: Nil if WomFloatType.isCoerceableFrom(ftype) => Success(WomIntegerType)
    case _ => Failure(new Exception(s"Unexpected 'ceil' arguments: $params (expects a single Float argument)"))
  }
  def basename(params: Seq[Try[WomType]]): Try[WomType] = params.toList match {
    case Success(fType) :: Nil if WomStringType.isCoerceableFrom(fType) => Success(WomStringType)
    case Success(fType) :: Success(sType) :: Nil if WomStringType.isCoerceableFrom(fType) && WomStringType.isCoerceableFrom(sType) => Success(WomStringType)
    case _ => Failure(new Exception(s"Unexpected basename arguments: $params"))
  }
  def transpose(params: Seq[Try[WomType]]): Try[WomType] = params.toList match {
    case Success(t @ WomArrayType(WomArrayType(_))) :: Nil => Success(t)
    case _ => Failure(new Exception(s"Unexpected transpose target: $params"))
  }
  def select_first(params: Seq[Try[WomType]]): Try[WomType] = extractSingleArgument("select_first", params) flatMap {
    case WomArrayType(WomOptionalType(innerType)) => Success(innerType)
    case other => Failure(new IllegalArgumentException(s"select_first failed. It expects an array of optional values but got ${other.toDisplayString}."))
  }
  def select_all(params: Seq[Try[WomType]]): Try[WomType] = extractSingleArgument("select_all", params) flatMap {
    case WomArrayType(WomOptionalType(innerType)) => Success(WomArrayType(innerType))
    case other => Failure(new IllegalArgumentException(s"select_all failed. It expects an array of optional values but got ${other.toDisplayString}."))
  }
  def defined(params: Seq[Try[WomType]]): Try[WomType] = extractSingleArgument("defined", params).map(_ => WomBooleanType)
  def zip(params: Seq[Try[WomType]]): Try[WomType] = {
    val badArgsFailure = Failure(new Exception(s"Unexpected zip parameters: $params"))
    WdlStandardLibraryFunctions.extractTwoParams(params, badArgsFailure) flatMap {
      case (arrayType1: WomArrayType, arrayType2: WomArrayType) => Success(WomArrayType(WomPairType(arrayType1.memberType, arrayType2.memberType)))
      case _ => badArgsFailure
    }
  }
  def cross(params: Seq[Try[WomType]]): Try[WomType] = zip(params)

  def flatten(params: Seq[Try[WomType]]): Try[WomType] = extractSingleArgument("flatten", params) flatMap {
    case WomArrayType(inner @ WomArrayType(_)) => Success(inner)
    case otherType => Failure(new Exception(s"flatten requires an Array[Array[_]] argument but instead got ${otherType.toDisplayString}"))
  }
}

case object NoFunctions extends WdlStandardLibraryFunctions {
  override def globHelper(pattern: String): Seq[String] = throw new NotImplementedError()
  override def readFile(path: String): String = throw new NotImplementedError()
  override def writeFile(path: String, content: String): Try[WomFile] = throw new NotImplementedError()
  override def stdout(params: Seq[Try[WomValue]]): Try[WomFile] = Failure(new NotImplementedError())
  override def stderr(params: Seq[Try[WomValue]]): Try[WomFile] = Failure(new NotImplementedError())
  override def read_json(params: Seq[Try[WomValue]]): Try[WomValue] = Failure(new NotImplementedError())
  override def write_tsv(params: Seq[Try[WomValue]]): Try[WomFile] = Failure(new NotImplementedError())
  override def write_json(params: Seq[Try[WomValue]]): Try[WomFile] = Failure(new NotImplementedError())
  override def size(params: Seq[Try[WomValue]]): Try[WomFloat] = Failure(new NotImplementedError())
  override def length(params: Seq[Try[WomValue]]): Try[WomInteger] = Failure(new NotImplementedError())
  override def flatten(params: Seq[Try[WomValue]]): Try[WomValue] = Failure(new NotImplementedError())
  override def sub(params: Seq[Try[WomValue]]): Try[WomString] = Failure(new NotImplementedError())
  override def range(params: Seq[Try[WomValue]]): Try[WomArray] = Failure(new NotImplementedError())
  override def transpose(params: Seq[Try[WomValue]]): Try[WomArray] = Failure(new NotImplementedError())
  override def select_first(params: Seq[Try[WomValue]]): Try[WomArray] = Failure(new NotImplementedError())
  override def select_all(params: Seq[Try[WomValue]]): Try[WomArray] = Failure(new NotImplementedError())
  override def zip(params: Seq[Try[WomValue]]): Try[WomArray] = Failure(new NotImplementedError())
  override def cross(params: Seq[Try[WomValue]]): Try[WomArray] = Failure(new NotImplementedError())
  override def basename(params: Seq[Try[WomValue]]): Try[WomString] = Failure(new NotImplementedError())
  override def floor(params: Seq[Try[WomValue]]): Try[WomInteger] = Failure(new NotImplementedError())
  override def round(params: Seq[Try[WomValue]]): Try[WomInteger] = Failure(new NotImplementedError())
  override def ceil(params: Seq[Try[WomValue]]): Try[WomInteger] = Failure(new NotImplementedError())
}
