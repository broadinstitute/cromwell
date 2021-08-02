package wdl.transforms.base.linking.expression.values

import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import common.validation.ErrorOr
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import wom.format.MemorySize
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{EvaluatedValue, ForCommandInstantiationOptions, ValueEvaluator}
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.shared.transforms.evaluation.values.EngineFunctions
import wdl4s.parser.MemoryUnit
import wom.expression.IoFunctionSet
import wom.types._
import wom.values.WomArray.WomArrayLike
import wom.values.{WomArray, WomBoolean, WomFloat, WomInteger, WomMap, WomObject, WomOptionalValue, WomPair, WomSingleFile, WomString, WomValue}
import wom.types.coercion.ops._
import wom.types.coercion.defaults._
import wom.types.coercion.WomTypeCoercer
import spray.json._
import wdl.model.draft3.elements.ExpressionElement
import wdl.shared.FileSizeLimitationConfig
import wdl.shared.model.expression.ValueEvaluation
import wom.CommandSetupSideEffectFile

import scala.concurrent.duration._
import scala.concurrent.Await
import scala.util.Try

object EngineFunctionEvaluators {
  private val fileSizeLimitationConfig = FileSizeLimitationConfig.fileSizeLimitationConfig

  implicit val stdoutFunctionEvaluator: ValueEvaluator[StdoutElement.type] = new ValueEvaluator[StdoutElement.type] {
    override def evaluateValue(a: StdoutElement.type,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomSingleFile]] =
      EvaluatedValue(WomSingleFile(ioFunctionSet.pathFunctions.stdout), Seq.empty).validNel
  }

  implicit val stderrFunctionEvaluator: ValueEvaluator[StderrElement.type] = new ValueEvaluator[StderrElement.type] {
    override def evaluateValue(a: StderrElement.type,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomSingleFile]] =
      EvaluatedValue(WomSingleFile(ioFunctionSet.pathFunctions.stderr), Seq.empty).validNel
  }

  private val ReadWaitTimeout = 300.seconds
  private def readFile(fileToRead: WomSingleFile, ioFunctionSet: IoFunctionSet, sizeLimit: Int) = {
    Try(Await.result(ioFunctionSet.readFile(fileToRead.value, Option(sizeLimit), failOnOverflow = true), ReadWaitTimeout))
  }

  implicit val readLinesFunctionEvaluator: ValueEvaluator[ReadLines] = new ValueEvaluator[ReadLines] {
    override def evaluateValue(a: ReadLines,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] = {
      processValidatedSingleValue[WomSingleFile, WomArray](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { fileToRead =>
        val tryResult = for {
          //validate
          read <- readFile(fileToRead, ioFunctionSet, fileSizeLimitationConfig.readLinesLimit)
          lines = read.split(System.lineSeparator)
        } yield EvaluatedValue(WomArray(lines map WomString.apply), Seq.empty)
        tryResult.toErrorOr.contextualizeErrors(s"""read_lines("${fileToRead.value}")""")
      }
    }
  }

  implicit val readTsvFunctionEvaluator: ValueEvaluator[ReadTsv] = new ValueEvaluator[ReadTsv] {
    override def evaluateValue(a: ReadTsv,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] = {
      processValidatedSingleValue[WomSingleFile, WomArray](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { fileToRead =>
        val tryResult = for {
          read <- readFile(fileToRead, ioFunctionSet, fileSizeLimitationConfig.readTsvLimit)
          tsv <- Try(WomArray.fromTsv(read))
        } yield EvaluatedValue(tsv, Seq.empty)
        tryResult.toErrorOr.contextualizeErrors(s"""read_tsv("${fileToRead.value}")""")
      }
    }
  }

  implicit val readMapFunctionEvaluator: ValueEvaluator[ReadMap] = new ValueEvaluator[ReadMap] {
    override def evaluateValue(a: ReadMap,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomMap]] = {
      processValidatedSingleValue[WomSingleFile, WomMap](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { fileToRead =>
        val tryResult = for {
          read <- readFile(fileToRead, ioFunctionSet, fileSizeLimitationConfig.readMapLimit)
          map <- WomMap.fromTsv(read)
        } yield EvaluatedValue(map, Seq.empty)
        tryResult.toErrorOr.contextualizeErrors(s"""read_map("${fileToRead.value}")""")
      }
    }
  }

  implicit val readObjectFunctionEvaluator: ValueEvaluator[ReadObject] = new ValueEvaluator[ReadObject] {
    override def evaluateValue(a: ReadObject,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomObject]] = {
      processValidatedSingleValue[WomSingleFile, WomObject](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { fileToRead =>
        val tryResult = for {
          read <- readFile(fileToRead, ioFunctionSet, fileSizeLimitationConfig.readObjectLimit)
          obj <- WomObject.fromTsv(read)
        } yield obj
        val rightSize: ErrorOr[WomObject] = tryResult.toErrorOr flatMap {
          case oneItem: Array[WomObject] if oneItem.length == 1 => oneItem.head.validNel
          case other: Array[WomObject] => s"Exactly 1 TSV object expected in input file (ie 2 lines: headers and data), but instead got an array of ${other.length} entries.".invalidNel
        }
        rightSize.map(EvaluatedValue(_, Seq.empty)).contextualizeErrors(s"""read_object("${fileToRead.value}")""")
      }
    }
  }

  implicit val readObjectsFunctionEvaluator: ValueEvaluator[ReadObjects] = new ValueEvaluator[ReadObjects] {
    override def evaluateValue(a: ReadObjects,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] = {
      processValidatedSingleValue[WomSingleFile, WomArray](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { fileToRead =>
        val tryResult = for {
          read <- readFile(fileToRead, ioFunctionSet, fileSizeLimitationConfig.readObjectLimit)
          objects <- WomObject.fromTsv(read)
        } yield WomArray(objects)

        tryResult.map(EvaluatedValue(_, Seq.empty)).toErrorOr.contextualizeErrors(s"""read_objects("${fileToRead.value}")""")
      }
    }
  }

  implicit val readJsonFunctionEvaluator: ValueEvaluator[ReadJson] = new ValueEvaluator[ReadJson] {
    override def evaluateValue(a: ReadJson,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomValue]] = {

      def convertJsonToWom(jsValue: JsValue): Try[WomValue] = {
        jsValue match {
          case _: JsNumber => WomIntegerType.coerceRawValue(jsValue).recoverWith { case _ => WomFloatType.coerceRawValue(jsValue) }
          case _: JsString => WomStringType.coerceRawValue(jsValue)
          case _: JsBoolean => WomBooleanType.coerceRawValue(jsValue)
          case _: JsArray => WomArrayType(WomAnyType).coerceRawValue(jsValue)
          case _ => WomObjectType.coerceRawValue(jsValue)
        }
      }

      def readJson(fileToRead: WomSingleFile): ErrorOr[EvaluatedValue[WomValue]] = {
        val tryResult: Try[WomValue] = for {
          read <- readFile(fileToRead, ioFunctionSet, fileSizeLimitationConfig.readJsonLimit)
          jsValue <- Try(read.parseJson)
          womValue <- convertJsonToWom(jsValue)
        } yield womValue

        tryResult.map(EvaluatedValue(_, Seq.empty)).toErrorOr.contextualizeErrors(s"""read_json("${fileToRead.value}")""")
      }

      def convertToSingleFile(womValue: WomValue): ErrorOr[WomSingleFile] = {
        if (womValue.coercionDefined[WomSingleFile]) womValue.coerceToType[WomSingleFile]
        else s"Expected File argument but got ${womValue.womType.stableName}".invalidNel
      }

      for {
        evaluatedValue <- a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        fileToRead <- convertToSingleFile(evaluatedValue.value)
        womValue <- readJson(fileToRead)
      } yield womValue
    }
  }

  implicit val readIntFunctionEvaluator: ValueEvaluator[ReadInt] = new ValueEvaluator[ReadInt] {
    override def evaluateValue(a: ReadInt,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomInteger]] = {
      processValidatedSingleValue[WomSingleFile, WomInteger](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { fileToRead =>
        val tryResult = for {
          read <- readFile(fileToRead, ioFunctionSet, fileSizeLimitationConfig.readIntLimit)
          asInt <- Try(read.trim.toInt)
        } yield WomInteger(asInt)
        tryResult.map(EvaluatedValue(_, Seq.empty)).toErrorOr.contextualizeErrors(s"""read_int("${fileToRead.value}")""")
      }
    }
  }

  implicit val readStringFunctionEvaluator: ValueEvaluator[ReadString] = new ValueEvaluator[ReadString] {
    override def evaluateValue(a: ReadString,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomString]] = {
      processValidatedSingleValue[WomSingleFile, WomString](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { fileToRead =>
        val tryResult = for {
          read <- readFile(fileToRead, ioFunctionSet, fileSizeLimitationConfig.readStringLimit)
        } yield WomString(read.trim)
        tryResult.map(EvaluatedValue(_, Seq.empty)).toErrorOr.contextualizeErrors(s"""read_string("${fileToRead.value}")""")
      }
    }
  }

  implicit val readFloatFunctionEvaluator: ValueEvaluator[ReadFloat] = new ValueEvaluator[ReadFloat] {
    override def evaluateValue(a: ReadFloat,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomFloat]] = {
      processValidatedSingleValue[WomSingleFile, WomFloat](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { fileToRead =>
        val tryResult = for {
          read <- readFile(fileToRead, ioFunctionSet, fileSizeLimitationConfig.readFloatLimit)
          asFloat <- Try(read.trim.toDouble)
        } yield WomFloat(asFloat)
        tryResult.map(EvaluatedValue(_, Seq.empty)).toErrorOr.contextualizeErrors(s"""read_float("${fileToRead.value}")""")
      }
    }
  }

  implicit val readBooleanFunctionEvaluator: ValueEvaluator[ReadBoolean] = new ValueEvaluator[ReadBoolean] {
    override def evaluateValue(a: ReadBoolean,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomBoolean]] = {
      processValidatedSingleValue[WomSingleFile, WomBoolean](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { fileToRead =>
        val tryResult = for {
          read <- readFile(fileToRead, ioFunctionSet, fileSizeLimitationConfig.readBoolLimit)
          asBool <- Try(read.trim.toBoolean)
        } yield WomBoolean(asBool)
        tryResult.map(EvaluatedValue(_, Seq.empty)).toErrorOr.contextualizeErrors(s"""read_boolean("${fileToRead.value}")""")
      }
    }
  }

  private val WriteWaitTimeout = 10.minutes
  private def writeContent(functionName: String, ioFunctionSet: IoFunctionSet, content: String): Try[WomSingleFile] = {
    import wom.values.HashableString
    Try(Await.result(ioFunctionSet.writeFile(s"${functionName}_${content.md5Sum}.tmp", content), WriteWaitTimeout))
  }

  implicit val writeLinesFunctionEvaluator: ValueEvaluator[WriteLines] = new ValueEvaluator[WriteLines] {
    override def evaluateValue(a: WriteLines,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomSingleFile]] = {
      val functionName = "write_lines"
      processValidatedSingleValue[WomArray, WomSingleFile](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { stringsToWrite =>
        val tryResult = for {
          serialized <- ValueEvaluation.serializeWomValue(functionName, stringsToWrite, defaultIfOptionalEmpty = WomArray(WomArrayType(WomStringType), Seq.empty))
          written <- writeContent(functionName, ioFunctionSet, serialized)
        } yield written

        tryResult.map(v => EvaluatedValue(v, Seq(CommandSetupSideEffectFile(v)))).toErrorOr.contextualizeErrors(s"""$functionName(...)""")
      } (coercer = WomArrayType(WomStringType))
    }
  }

  implicit val writeTsvFunctionEvaluator: ValueEvaluator[WriteTsv] = new ValueEvaluator[WriteTsv] {
    override def evaluateValue(a: WriteTsv,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomSingleFile]] = {
      val functionName = "write_tsv"
      processValidatedSingleValue[WomArray, WomSingleFile](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { objectToWrite =>

        val tryResult = for {
          serialized <- ValueEvaluation.serializeWomValue(functionName, objectToWrite, defaultIfOptionalEmpty = WomArray(WomArrayType(WomStringType), List.empty[WomValue]))
          written <- writeContent(functionName, ioFunctionSet, serialized)
        } yield written

        tryResult.map(v => EvaluatedValue(v, Seq(CommandSetupSideEffectFile(v)))).toErrorOr.contextualizeErrors(s"""$functionName(...)""")
      } (coercer = WomArrayType(WomAnyType))
    }
  }

  implicit val writeMapFunctionEvaluator: ValueEvaluator[WriteMap] = new ValueEvaluator[WriteMap] {
    override def evaluateValue(a: WriteMap,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomSingleFile]] = {
      val functionName = "write_map"
      processValidatedSingleValue[WomMap, WomSingleFile](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { mapToWrite: WomMap =>
        val tryResult = for {
          serialized <- ValueEvaluation.serializeWomValue(functionName, mapToWrite, defaultIfOptionalEmpty = WomMap(Map.empty))
          written <- writeContent(functionName, ioFunctionSet, serialized)
        } yield written

        tryResult.map(v => EvaluatedValue(v, Seq(CommandSetupSideEffectFile(v)))).toErrorOr.contextualizeErrors(s"""$functionName(...)""")
      }
    }
  }

  implicit val writeObjectFunctionEvaluator: ValueEvaluator[WriteObject] = new ValueEvaluator[WriteObject] {
    override def evaluateValue(a: WriteObject,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomSingleFile]] = {
      val functionName = "write_object"
      processValidatedSingleValue[WomObject, WomSingleFile](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { objectToWrite: WomObject =>
        val tryResult = for {
          serialized <- ValueEvaluation.serializeWomValue(functionName, objectToWrite, defaultIfOptionalEmpty = WomObject(Map.empty))
          written <- writeContent(functionName, ioFunctionSet, serialized)
        } yield written

        tryResult.map(v => EvaluatedValue(v, Seq(CommandSetupSideEffectFile(v)))).toErrorOr.contextualizeErrors(s"""$functionName(...)""")
      }
    }
  }

  implicit val writeObjectsFunctionEvaluator: ValueEvaluator[WriteObjects] = new ValueEvaluator[WriteObjects] {
    override def evaluateValue(a: WriteObjects,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomSingleFile]] = {
      val functionName = "write_objects"
      processValidatedSingleValue[WomArray, WomSingleFile](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { objectToWrite =>
        val tryResult = for {
          serialized <- ValueEvaluation.serializeWomValue(functionName, objectToWrite, defaultIfOptionalEmpty = WomArray(List(WomObject(Map.empty))))
          written <- writeContent(functionName, ioFunctionSet, serialized)
        } yield written

        tryResult.map(v => EvaluatedValue(v, Seq(CommandSetupSideEffectFile(v)))).toErrorOr.contextualizeErrors(s"""$functionName(...)""")
      } (coercer = WomArrayType(WomObjectType))
    }
  }

  implicit val writeJsonFunctionEvaluator: ValueEvaluator[WriteJson] = new ValueEvaluator[WriteJson] {
    override def evaluateValue(a: WriteJson,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomSingleFile]] = {
      val functionName = "write_json"

      def convertToSingleFile(objectToWrite: WomValue): ErrorOr[EvaluatedValue[WomSingleFile]] = {
        val serialized = ValueEvaluation.valueToJson(objectToWrite)
        val tryResult = for {
          written <- writeContent(functionName, ioFunctionSet, serialized.compactPrint)
        } yield written

        tryResult.map(v => EvaluatedValue(v, Seq(CommandSetupSideEffectFile(v)))).toErrorOr.contextualizeErrors(s"""$functionName(...)""")
      }

      def evaluateParam(womValue: WomValue): ErrorOr[EvaluatedValue[WomSingleFile]] = {
        womValue match {
          case WomBoolean(_) | WomString(_) | WomInteger(_) | WomFloat(_) | WomPair(_, _) => convertToSingleFile(womValue)
          case v if v.coercionDefined[WomObject] => v.coerceToType[WomObject].flatMap(convertToSingleFile)
          case v if v.coercionDefined[WomArray] => v.coerceToType[WomArray].flatMap(convertToSingleFile)
          case _ => (s"The '$functionName' method expects one of 'Boolean', 'String', 'Integer', 'Float', 'Object', 'Pair[_, _]', " +
            s"'Map[_, _] or 'Array[_]' argument but instead got '${womValue.womType.friendlyName}'.").invalidNel
        }
      }

      val evaluatedSingleFile: ErrorOr[(EvaluatedValue[WomSingleFile], Seq[CommandSetupSideEffectFile])] = for {
        evaluatedValue <- a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
        evaluatedSingleFile <- evaluateParam(evaluatedValue.value)
      } yield (evaluatedSingleFile, evaluatedValue.sideEffectFiles)

      evaluatedSingleFile map {
        case (result, previousSideEffectFiles) => result.copy(sideEffectFiles = result.sideEffectFiles ++ previousSideEffectFiles)
      }
    }
  }

  implicit val rangeFunctionEvaluator: ValueEvaluator[Range] = new ValueEvaluator[Range] {
    override def evaluateValue(a: Range,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] = {
      processValidatedSingleValue[WomInteger, WomArray](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { integer =>
        val array = WomArray(
          womType = WomArrayType(WomIntegerType, guaranteedNonEmpty = integer.value > 0),
          value = (0 until integer.value).map(WomInteger)
        )
        EvaluatedValue(array, Seq.empty).validNel
      }
    }
  }

  implicit val transposeFunctionEvaluator: ValueEvaluator[Transpose] = new ValueEvaluator[Transpose] {
    override def evaluateValue(a: Transpose,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] = {
      processValidatedSingleValue[WomArray, WomArray](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { array =>
        EngineFunctions.transpose(array).map(EvaluatedValue(_, Seq.empty)).toErrorOr
      }
    }
  }

  implicit val lengthFunctionEvaluator: ValueEvaluator[Length] = new ValueEvaluator[Length] {
    override def evaluateValue(a: Length,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomInteger]] = {
      processValidatedSingleValue[WomArray, WomInteger](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { a =>
        EvaluatedValue(WomInteger(a.value.size), Seq.empty).validNel
      }
    }
  }

  implicit val flattenFunctionEvaluator: ValueEvaluator[Flatten] = new ValueEvaluator[Flatten] {
    override def evaluateValue(a: Flatten,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] = {
      def flatValues(v: WomValue): ErrorOr[Seq[WomValue]] = v match {
        case WomArrayLike(arrayLike) => arrayLike.value.validNel
        case other => s"inner item ${other.toWomString} was not an array-like".invalidNel
      }

      processValidatedSingleValue[WomArray, WomArray](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { array =>
        val expandedValidation = array.value.toList.traverse{ flatValues }
        expandedValidation map { expanded => EvaluatedValue(WomArray(expanded.flatten), Seq.empty) }
      } (coercer = WomArrayType(WomArrayType(WomAnyType)))
    }
  }

  implicit val selectFirstFunctionEvaluator: ValueEvaluator[SelectFirst] = new ValueEvaluator[SelectFirst] {
    override def evaluateValue(a: SelectFirst,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomValue]] = {
      processValidatedSingleValue[WomArray, WomValue](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { array =>
        val firstValue = array.value collectFirst {
          case WomOptionalValue(_, Some(yay)) => yay
        }

        firstValue match {
          case Some(first) => EvaluatedValue(first, Seq.empty).validNel
          case None => s"select_first was called with ${array.size} empty values. We needed at least one to be filled.".invalidNel
        }
      } (coercer = WomArrayType(WomOptionalType(WomAnyType)))
    }
  }

  implicit val selectAllFunctionEvaluator: ValueEvaluator[SelectAll] = new ValueEvaluator[SelectAll] {
    override def evaluateValue(a: SelectAll,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] = {
      processValidatedSingleValue[WomArray, WomArray](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { array =>
        val goodValues = array.value.collect {
          case WomOptionalValue.Flattened(Some(value)) => value
        }
        EvaluatedValue(WomArray(goodValues), Seq.empty).validNel
      } (coercer = WomArrayType(WomOptionalType(WomAnyType)))
    }
  }

  implicit val definedFunctionEvaluator: ValueEvaluator[Defined] = new ValueEvaluator[Defined] {
    override def evaluateValue(a: Defined,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomBoolean]] = {
      processValidatedSingleValue[WomOptionalValue, WomBoolean](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { optionalValue =>
        EvaluatedValue(WomBoolean(optionalValue.value.isDefined), Seq.empty).validNel
      }
    }
  }

  implicit val floorFunctionEvaluator: ValueEvaluator[Floor] = new ValueEvaluator[Floor] {
    override def evaluateValue(a: Floor,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomInteger]] = {
      processValidatedSingleValue[WomFloat, WomInteger](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { float =>
        EvaluatedValue(WomInteger(math.floor(float.value).toInt), Seq.empty).validNel
      }
    }
  }

  implicit val ceilFunctionEvaluator: ValueEvaluator[Ceil] = new ValueEvaluator[Ceil] {
    override def evaluateValue(a: Ceil,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomInteger]] = {
      processValidatedSingleValue[WomFloat, WomInteger](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { float =>
        EvaluatedValue(WomInteger(math.ceil(float.value).toInt), Seq.empty).validNel
      }
    }
  }

  implicit val roundFunctionEvaluator: ValueEvaluator[Round] = new ValueEvaluator[Round] {
    override def evaluateValue(a: Round,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomInteger]] = {
      processValidatedSingleValue[WomFloat, WomInteger](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { float =>
        EvaluatedValue(WomInteger(math.round(float.value).toInt), Seq.empty).validNel
      }
    }
  }

  implicit val globFunctionValueEvaluator: ValueEvaluator[Glob] = new ValueEvaluator[Glob] {
    /**
      * Evaluate a value from an A
      *
      * @param a                              The A to evaluate
      * @param inputs                         Evaluation inputs
      * @param ioFunctionSet                  IO functions to use
      * @param forCommandInstantiationOptions Supplied only if we're evaluating this A as part of command instantiation.
      * @return An evaluated value set - the value itself and any files which were produced as part of the evaluation.
      */
    override def evaluateValue(a: Glob, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
      processValidatedSingleValue[WomString, WomArray](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { globString =>
        for {
          globbed <- Try(Await.result(ioFunctionSet.glob(globString.valueString), ReadWaitTimeout)).toErrorOr
          files = globbed map WomSingleFile
          array = WomArray(files)
        } yield EvaluatedValue(array, Seq.empty)
      }
    }
  }

  implicit val sizeFunctionEvaluator: ValueEvaluator[Size] = new ValueEvaluator[Size] {
    override def evaluateValue(a: Size,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomFloat]] = {
      // Inner function: get the memory unit from the second (optional) parameter
      def toUnit(womValue: WomValue): ErrorOr[MemoryUnit] = Try(MemoryUnit.fromSuffix(womValue.valueString)).toErrorOr

      // Inner function: is this a file type, or an optional containing a file type?
      def isOptionalOfFileType(womType: WomType): Boolean = womType match {
        case f if WomSingleFileType.isCoerceableFrom(f) => true
        case WomOptionalType(inner) => isOptionalOfFileType(inner)
        case _ => false
      }

      // Inner function: Get the file size, allowing for unpacking of optionals and arrays
      def optionalSafeFileSize(value: WomValue): ErrorOr[Long] = value match {
        case f if f.isInstanceOf[WomSingleFile] || WomSingleFileType.isCoerceableFrom(f.womType) =>
          f.coerceToType[WomSingleFile] flatMap { file => Try(Await.result(ioFunctionSet.size(file.valueString), Duration.Inf)).toErrorOr }
        case WomOptionalValue(f, Some(o)) if isOptionalOfFileType(f) => optionalSafeFileSize(o)
        case WomOptionalValue(f, None) if isOptionalOfFileType(f) => 0l.validNel
        case WomArray(WomArrayType(womType), values) if isOptionalOfFileType(womType) => values.toList.traverse(optionalSafeFileSize).map(_.sum)
        case _ => s"The 'size' method expects a 'File', 'File?', 'Array[File]' or Array[File?] argument but instead got ${value.womType.stableName}.".invalidNel
      }

      // Inner function: get the file size and convert into the requested memory unit
      def fileSize(womValue: ErrorOr[EvaluatedValue[_ <: WomValue]], convertToOption: Option[ErrorOr[EvaluatedValue[_ <: WomValue]]]): ErrorOr[EvaluatedValue[WomFloat]] = {
        val convertTo: ErrorOr[EvaluatedValue[_ <: WomValue]] = convertToOption.getOrElse(EvaluatedValue(WomString("B"), Seq.empty).validNel)
        for {
          value <- womValue
          evaluatedUnitValue <- convertTo
          convertToUnit <- toUnit(evaluatedUnitValue.value)
          fileSize <- optionalSafeFileSize(value.value)
        } yield EvaluatedValue(WomFloat(MemorySize(fileSize.toDouble, MemoryUnit.Bytes).to(convertToUnit).amount), value.sideEffectFiles ++ evaluatedUnitValue.sideEffectFiles)
      }

      val evaluatedFileValidation: ErrorOr[EvaluatedValue[_ <: WomValue]] = a.file.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)
      fileSize(evaluatedFileValidation, a.unit map (_.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)))
    }
  }

  implicit val basenameFunctionEvaluator: ValueEvaluator[Basename] = new ValueEvaluator[Basename] {
    override def evaluateValue(a: Basename,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomString]] = {
      def simpleBasename(fileNameAsString: WomString) = fileNameAsString.valueString.split('/').last

      a.suffixToRemove match {
        case None => processValidatedSingleValue[WomString, WomString](a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { str =>
          EvaluatedValue(WomString(simpleBasename(str)), Seq.empty).validNel
        }
        case Some(suffixToRemove) => processTwoValidatedValues[WomString, WomString, WomString](
          a.param.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions),
          suffixToRemove.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { (name, suffix) =>
            EvaluatedValue(WomString(simpleBasename(name).stripSuffix(suffix.valueString)), Seq.empty).validNel
          }
      }
    }
  }

  implicit val zipFunctionEvaluator: ValueEvaluator[Zip] = new ValueEvaluator[Zip] {
    override def evaluateValue(a: Zip,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] = {
      processTwoValidatedValues[WomArray, WomArray, WomArray](a.arg1.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions), a.arg2.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { (arr1, arr2) =>
        if (arr1.size == arr2.size) {
          val pairs = arr1.value.zip(arr2.value) map { case (a, b) => WomPair(a, b) }
          EvaluatedValue(WomArray(WomArrayType(WomPairType(arr1.arrayType.memberType, arr2.arrayType.memberType)), pairs), Seq.empty).validNel
        } else {
          s"Mismatching array sizes for zip function: ${arr1.size} vs ${arr2.size}".invalidNel
        }
      }
    }
  }

  implicit val crossFunctionEvaluator: ValueEvaluator[Cross] = new ValueEvaluator[Cross] {
    override def evaluateValue(a: Cross,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] = {
      processTwoValidatedValues[WomArray, WomArray, WomArray](a.arg1.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions), a.arg2.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { (arr1, arr2) =>
        val pairs = for {
          a <- arr1.value
          b <- arr2.value
        } yield WomPair(a, b)
        EvaluatedValue(WomArray(WomArrayType(WomPairType(arr1.arrayType.memberType, arr2.arrayType.memberType)), pairs), Seq.empty).validNel
      }
    }
  }

  implicit val prefixFunctionEvaluator: ValueEvaluator[Prefix] = new ValueEvaluator[Prefix] {
    override def evaluateValue(a: Prefix,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomArray]] = {
      processTwoValidatedValues[WomString, WomArray, WomArray](a.prefix.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions), a.array.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { (prefix, array) =>
        EvaluatedValue(WomArray(array.value.map(value => WomString(prefix.value + value.valueString))), Seq.empty).validNel
      }
    }
  }

  implicit val subFunctionEvaluator: ValueEvaluator[Sub] = new ValueEvaluator[Sub] {
    override def evaluateValue(a: Sub,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[WomString]] = {
      processThreeValidatedValues[WomString, WomString, WomString, WomString](
        a.input.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions),
        a.pattern.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions),
        a.replace.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) { (input, pattern, replace) =>
        ErrorOr(EvaluatedValue(WomString(pattern.valueString.r.replaceAllIn(input.valueString, replace.valueString)), Seq.empty))
      }
    }
  }

  def processValidatedSingleValue[A <: WomValue, B <: WomValue](arg: ErrorOr[EvaluatedValue[_]])
                                                                       (f: A => ErrorOr[EvaluatedValue[B]])
                                                                       (implicit coercer: WomTypeCoercer[A]): ErrorOr[EvaluatedValue[B]] = {
    arg flatMap {
      case EvaluatedValue(a: WomValue, previousSideEffectFiles) if a.coercionDefined[A] => a.coerceToType[A] flatMap { f.apply } map { result => result.copy(sideEffectFiles = result.sideEffectFiles ++ previousSideEffectFiles) }
      case other => s"Expected ${coercer.toDisplayString} argument but got ${other.value.womType.stableName}".invalidNel
    }
  }

  def processTwoValidatedValues[A <: WomValue, B <: WomValue, R <: WomValue](arg1: ErrorOr[EvaluatedValue[_ <: WomValue]], arg2: ErrorOr[EvaluatedValue[_ <: WomValue]])
                                                                     (f: (A, B) => ErrorOr[EvaluatedValue[R]])
                                                                     (implicit coercerA: WomTypeCoercer[A],
                                                                      coercerB: WomTypeCoercer[B]): ErrorOr[EvaluatedValue[R]] = {
    (arg1, arg2) flatMapN {
      case (EvaluatedValue(a: WomValue, previousSideEffectFilesA), EvaluatedValue(b: WomValue, previousSideEffectFilesB)) if a.coercionDefined[A] && b.coercionDefined[B] =>
        (a.coerceToType[A], b.coerceToType[B]) flatMapN { f.apply } map { result => result.copy(sideEffectFiles = result.sideEffectFiles ++ previousSideEffectFilesA ++ previousSideEffectFilesB) }
      case (otherA, otherB) => s"Expected (${coercerA.toDisplayString}, ${coercerB.toDisplayString}) argument but got (${otherA.value.womType.stableName}, ${otherB.value.womType.stableName})".invalidNel
    }
  }

  private def processThreeValidatedValues[A <: WomValue, B <: WomValue, C <: WomValue, R <: WomValue](arg1: ErrorOr[EvaluatedValue[_ <: WomValue]], arg2: ErrorOr[EvaluatedValue[_ <: WomValue]], arg3: ErrorOr[EvaluatedValue[_ <: WomValue]])
                                                                                      (f: (A, B, C) => ErrorOr[EvaluatedValue[R]])
                                                                                      (implicit coercerA: WomTypeCoercer[A],
                                                                                       coercerB: WomTypeCoercer[B],
                                                                                       coercerC: WomTypeCoercer[C]): ErrorOr[EvaluatedValue[R]] = {
    (arg1, arg2, arg3) flatMapN {
      case (EvaluatedValue(a, previousSideEffectFilesA), EvaluatedValue(b, previousSideEffectFilesB), EvaluatedValue(c, previousSideEffectFilesC)) if a.coercionDefined[A] && b.coercionDefined[B] && c.coercionDefined[C] =>
        (a.coerceToType[A], b.coerceToType[B], c.coerceToType[C]) flatMapN { f.apply } map { result => result.copy(sideEffectFiles = result.sideEffectFiles ++ previousSideEffectFilesA ++ previousSideEffectFilesB ++ previousSideEffectFilesC) }
      case (otherA, otherB, otherC) => s"Expected (${coercerA.toDisplayString}, ${coercerB.toDisplayString}, ${coercerB.toDisplayString}) argument but got (${otherA.value.womType.stableName}, ${otherB.value.womType.stableName}, ${otherC.value.womType.stableName})".invalidNel
    }
  }
}
