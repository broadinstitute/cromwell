package wdl.transforms.v1_2.linking.expression.values

import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.google.re2j.{Pattern => RE2JPattern}
import common.validation.ErrorOr._
import common.validation.Validation._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{
  Equals,
  Keys,
  Length,
  PrimitiveLiteralExpressionElement,
  Size,
  SubPosix,
  WriteTsv
}
import wdl.model.draft3.graph.expression.{EvaluatedValue, ForCommandInstantiationOptions, ValueEvaluator}
import wdl.model.v1_2.elements.V1_2ExpressionElement._
import wdl.shared.FileSizeLimitationConfig
import wdl.shared.model.expression.ValueEvaluation
import wdl4s.parser.MemoryUnit
import wom.CommandSetupSideEffectFile
import wom.expression.IoFunctionSet
import wom.format.MemorySize
import wom.types._
import wom.values._
import wom.values.WomArray.WomArrayLike

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.util.Try

/**
 * Value evaluation for the operators and functions added in WDL 1.2, and for the WDL 1.0/1.1 functions whose accepted
 * argument types or semantics changed in WDL 1.2.
 */
object V1_2FunctionValueEvaluators {

  private val IoTimeout = 10.minutes
  private val fileSizeLimitationConfig = FileSizeLimitationConfig.fileSizeLimitationConfig

  def evaluateValue(a: V1_2FunctionCallElement,
                    inputs: Map[String, WomValue],
                    ioFunctionSet: IoFunctionSet,
                    forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
  )(implicit expressionEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
    def evaluate(e: ExpressionElement): ErrorOr[EvaluatedValue[_ <: WomValue]] =
      expressionEvaluator.evaluateValue(e, inputs, ioFunctionSet, forCommandInstantiationOptions)(expressionEvaluator)

    // Evaluate all parameters, then compute the result from their values. Side effect files are carried along.
    def withParams(f: PartialFunction[List[WomValue], ErrorOr[WomValue]]): ErrorOr[EvaluatedValue[_ <: WomValue]] =
      a.params.toList.traverse(evaluate) flatMap { evaluated =>
        val values = evaluated.map(_.value)
        val sideEffectFiles = evaluated.flatMap(_.sideEffectFiles)
        f.applyOrElse(values,
                      (other: List[WomValue]) =>
                        s"Invalid arguments: (${other.map(_.womType.stableName).mkString(", ")})".invalidNel
        ).map(EvaluatedValue(_, sideEffectFiles))
      }

    a match {
      case _: Exponentiation =>
        withParams { case List(left, right) => exponentiation(left, right) }

      case _: Find =>
        withParams { case List(input, pattern) =>
          compile(pattern.valueString) map { regex =>
            val matcher = regex.matcher(input.valueString)
            if (matcher.find()) WomOptionalValue(WomString(matcher.group()))
            else WomOptionalValue(WomStringType, None)
          }
        }

      case _: Matches =>
        withParams { case List(input, pattern) =>
          compile(pattern.valueString).map(regex => WomBoolean(regex.matcher(input.valueString).find()))
        }

      case _: JoinPaths =>
        withParams {
          case List(WomArrayLike(paths)) => joinPaths(paths.value.map(unwrap(_).valueString).toList)
          case List(base, WomArrayLike(relative)) =>
            joinPaths(unwrap(base).valueString :: relative.value.map(unwrap(_).valueString).toList)
          case List(base, relative) => joinPaths(List(unwrap(base).valueString, unwrap(relative).valueString))
        }

      case ReadTsvWithHeader(_, header, fieldNames) =>
        withParams { case file :: WomBoolean(hasHeader) :: names =>
          val explicitNames = names.headOption.map {
            case WomArrayLike(array) => array.value.map(_.valueString).toList
            case other => List(other.valueString)
          }
          readTsv(unwrap(file).valueString, hasHeader, explicitNames, ioFunctionSet, header, fieldNames.isDefined)
        }

      case _: Contains =>
        withParams { case List(WomArrayLike(array), value) =>
          WomBoolean(array.value.exists(element => valuesEqual(element, value))).validNel
        }

      case _: Chunk =>
        withParams {
          case List(WomArrayLike(array), WomInteger(size)) if size > 0 =>
            val chunks = array.value.grouped(size).toList.map(chunk => WomArray(array.arrayType, chunk))
            WomArray(WomArrayType(array.arrayType), chunks).validNel
          case List(_, WomInteger(size)) => s"chunk size must be greater than 0 but was $size".invalidNel
        }

      case SelectFirstWithDefault(array, default) =>
        evaluate(array) flatMap { evaluatedArray =>
          val firstDefined = evaluatedArray.value match {
            case WomArrayLike(womArray) =>
              womArray.value.collectFirst {
                case WomOptionalValue.Flattened(Some(value)) => value
                case value if !value.isInstanceOf[WomOptionalValue] => value
              }
            case _ => None
          }
          firstDefined match {
            case Some(value) => EvaluatedValue(value, evaluatedArray.sideEffectFiles).validNel
            case None =>
              evaluate(default).map(d => EvaluatedValue(d.value, evaluatedArray.sideEffectFiles ++ d.sideEffectFiles))
          }
        }

      case _: ContainsKey =>
        withParams { case List(collection, key) => WomBoolean(containsKey(collection, key)).validNel }

      case _: Values =>
        withParams { case List(WomMap(WomMapType(_, valueType), map)) =>
          WomArray(WomArrayType(valueType), map.values.toList).validNel
        }

      case WriteTsvWithHeader(_, _, columnNames) =>
        withParams { case array :: WomBoolean(header) :: names =>
          val explicitNames = names.headOption.map {
            case WomArrayLike(namesArray) => namesArray.value.map(_.valueString).toList
            case other => List(other.valueString)
          }
          if (!header && columnNames.isDefined)
            "write_tsv: column names can only be specified when the header argument is true".invalidNel
          else writeTsv(array, header, explicitNames, ioFunctionSet)
        } map { evaluated =>
          evaluated.copy(sideEffectFiles = evaluated.sideEffectFiles ++ sideEffectFile(evaluated.value))
        }

      case _: EnvironmentVariableValue =>
        withParams { case List(value) => WomString(shellQuote(serializeForEnvironment(value))).validNel }
    }
  }

  // WDL 1.2: `write_tsv` also accepts an array of structs
  def evaluateWriteTsv(a: WriteTsv,
                       inputs: Map[String, WomValue],
                       ioFunctionSet: IoFunctionSet,
                       forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
  )(implicit expressionEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] =
    expressionEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(
      expressionEvaluator
    ) flatMap { evaluated =>
      writeTsv(evaluated.value, header = false, None, ioFunctionSet) map { file =>
        EvaluatedValue(file, evaluated.sideEffectFiles ++ sideEffectFile(file))
      }
    }

  /**
   * WDL 1.2: either or both operands of `==` and `!=` may be optional, and `None` is equal only to itself. When neither
   * operand is `None`, the WDL 1.0/1.1 equality rules (`fallback`) apply.
   */
  def evaluateEquality(left: ExpressionElement,
                       right: ExpressionElement,
                       negate: Boolean,
                       inputs: Map[String, WomValue],
                       ioFunctionSet: IoFunctionSet,
                       forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
  )(fallback: => ErrorOr[EvaluatedValue[_ <: WomValue]])(implicit
    expressionEvaluator: ValueEvaluator[ExpressionElement]
  ): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
    def evaluate(e: ExpressionElement) =
      expressionEvaluator.evaluateValue(e, inputs, ioFunctionSet, forCommandInstantiationOptions)(expressionEvaluator)
    def isNone(value: WomValue): Boolean = value match {
      case WomOptionalValue(_, None) => true
      case WomOptionalValue(_, Some(inner)) => isNone(inner)
      case _ => false
    }

    (evaluate(left), evaluate(right)) flatMapN { (l, r) =>
      val sideEffectFiles = l.sideEffectFiles ++ r.sideEffectFiles
      if (isNone(l.value) || isNone(r.value)) {
        val equal = isNone(l.value) && isNone(r.value)
        EvaluatedValue(WomBoolean(equal != negate), sideEffectFiles).validNel
      } else
        (unwrap(l.value), unwrap(r.value)) match {
          // Compare the defined values with the WDL 1.0/1.1 rules (e.g. Int and Float are compared numerically)
          case (lp: WomPrimitive, rp: WomPrimitive) if l.value != lp || r.value != rp =>
            evaluate(Equals(PrimitiveLiteralExpressionElement(lp), PrimitiveLiteralExpressionElement(rp))) map {
              case EvaluatedValue(WomBoolean(equal), _) => EvaluatedValue(WomBoolean(equal != negate), sideEffectFiles)
              case other => other
            }
          case _ => fallback
        }
    }
  }

  // WDL 1.2: `keys` also accepts a Struct or an Object
  def evaluateKeys(a: Keys,
                   inputs: Map[String, WomValue],
                   ioFunctionSet: IoFunctionSet,
                   forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
  )(implicit expressionEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] =
    expressionEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(
      expressionEvaluator
    ) flatMap { evaluated =>
      unwrap(evaluated.value) match {
        case WomMap(WomMapType(keyType, _), values) =>
          EvaluatedValue(WomArray(WomArrayType(keyType), values.keys.toList), evaluated.sideEffectFiles).validNel
        case WomObject(values, WomCompositeType(typeMap, _)) =>
          // Struct keys are returned in the order of the struct definition
          val orderedKeys = typeMap.keys.toList.filter(values.contains)
          EvaluatedValue(WomArray(WomArrayType(WomStringType), orderedKeys.map(WomString)),
                         evaluated.sideEffectFiles
          ).validNel
        case WomObject(values, _) =>
          EvaluatedValue(WomArray(WomArrayType(WomStringType), values.keys.toList.map(WomString)),
                         evaluated.sideEffectFiles
          ).validNel
        case other =>
          s"Invalid call of 'keys' on parameter of type '${other.womType.stableName}' (expected a map, struct or object)".invalidNel
      }
    }

  // WDL 1.2: `length` also accepts a Map, an Object or a String
  def evaluateLength(a: Length,
                     inputs: Map[String, WomValue],
                     ioFunctionSet: IoFunctionSet,
                     forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
  )(implicit expressionEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] =
    expressionEvaluator.evaluateValue(a.param, inputs, ioFunctionSet, forCommandInstantiationOptions)(
      expressionEvaluator
    ) flatMap { evaluated =>
      val length: ErrorOr[Int] = unwrap(evaluated.value) match {
        case WomArrayLike(array) => array.value.size.validNel
        case WomMap(_, map) => map.size.validNel
        case WomObject(values, _) => values.size.validNel
        case WomString(s) => s.codePointCount(0, s.length).validNel
        case other =>
          s"Invalid call of 'length' on parameter of type '${other.womType.stableName}' (expected an array, map, object or string)".invalidNel
      }
      length.map(l => EvaluatedValue(WomInteger(l), evaluated.sideEffectFiles))
    }

  // WDL 1.2: `size` also accepts any compound value containing files
  def evaluateSize(a: Size,
                   inputs: Map[String, WomValue],
                   ioFunctionSet: IoFunctionSet,
                   forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
  )(implicit expressionEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
    def evaluate(e: ExpressionElement) =
      expressionEvaluator.evaluateValue(e, inputs, ioFunctionSet, forCommandInstantiationOptions)(expressionEvaluator)

    val unitValidation: ErrorOr[Option[EvaluatedValue[_ <: WomValue]]] = a.unit.traverse(evaluate)
    (evaluate(a.file), unitValidation) flatMapN { (value, unitValue) =>
      for {
        unit <- Try(MemoryUnit.fromSuffix(unitValue.map(_.value.valueString).getOrElse("B"))).toErrorOr
        files <- filesToSize(value.value)
        bytes <- Try(Await.result(ioFunctionSet.parallelSize(files.map(_.valueString)), IoTimeout)).toErrorOr
      } yield EvaluatedValue(
        WomFloat(MemorySize(bytes.toDouble, MemoryUnit.Bytes).to(unit).amount),
        value.sideEffectFiles ++ unitValue.toSeq.flatMap(_.sideEffectFiles)
      )
    }
  }

  // WDL 1.2: `sub` supports `\1` through `\9` as references to capturing groups in the replacement string
  def evaluateSub(a: SubPosix,
                  inputs: Map[String, WomValue],
                  ioFunctionSet: IoFunctionSet,
                  forCommandInstantiationOptions: Option[ForCommandInstantiationOptions]
  )(implicit expressionEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] =
    List(a.input, a.pattern, a.replace)
      .traverse(
        expressionEvaluator.evaluateValue(_, inputs, ioFunctionSet, forCommandInstantiationOptions)(expressionEvaluator)
      )
      .flatMap {
        case evaluated @ List(input, pattern, replace) =>
          compile(pattern.value.valueString) map { regex =>
            val replaced =
              regex.matcher(input.value.valueString).replaceAll(toRe2jReplacement(replace.value.valueString))
            EvaluatedValue(WomString(replaced), evaluated.flatMap(_.sideEffectFiles))
          }
        case _ => "sub expects exactly 3 arguments".invalidNel
      }

  /**
   * Converts a WDL replacement string (where `\1`..`\9` are group references and other characters are literal) into an
   * RE2J replacement string (where `$1`..`$9` are group references and `\` escapes the next character).
   */
  def toRe2jReplacement(wdlReplacement: String): String = {
    val builder = new StringBuilder
    var i = 0
    while (i < wdlReplacement.length) {
      val c = wdlReplacement.charAt(i)
      if (
        c == '\\' && i + 1 < wdlReplacement.length && wdlReplacement.charAt(i + 1).isDigit &&
        wdlReplacement.charAt(i + 1) != '0'
      ) {
        builder.append('$').append(wdlReplacement.charAt(i + 1))
        i += 2
      } else {
        if (c == '\\' || c == '$') builder.append('\\')
        builder.append(c)
        i += 1
      }
    }
    builder.toString
  }

  // The spec requires POSIX ERE. RE2 is used, as for `sub` in WDL 1.1, since it supports the ERE syntax.
  private def compile(pattern: String): ErrorOr[RE2JPattern] =
    Try(RE2JPattern.compile(pattern)).toErrorOr

  private def unwrap(value: WomValue): WomValue = value match {
    case WomOptionalValue(_, Some(inner)) => unwrap(inner)
    case other => other
  }

  private def exponentiation(left: WomValue, right: WomValue): ErrorOr[WomValue] = (unwrap(left), unwrap(right)) match {
    case (WomInteger(base), WomInteger(exponent)) if exponent >= 0 =>
      val result = BigInt(base).pow(exponent)
      if (result.isValidInt) WomInteger(result.toInt).validNel
      else s"Integer overflow evaluating $base ** $exponent".invalidNel
    case (WomInteger(base), WomInteger(exponent)) =>
      s"Cannot evaluate $base ** $exponent: an Int cannot be raised to a negative Int power. Use a Float instead".invalidNel
    case (l, r) =>
      (toDouble(l), toDouble(r)) mapN { (base, exponent) => WomFloat(math.pow(base, exponent)) }
  }

  private def toDouble(value: WomValue): ErrorOr[Double] = value match {
    case WomInteger(i) => i.toDouble.validNel
    case WomFloat(f) => f.validNel
    case other => s"Expected a numeric value but got ${other.womType.stableName}".invalidNel
  }

  /**
   * TODO: The spec defines `join_paths` in terms of the execution environment's filesystem: paths are joined with `/`
   *  and the result is not checked for existence. Cloud URLs (e.g. `gs://bucket/dir`) are joined the same way.
   *  Suggestion: once Directory values are canonicalized by the backends (WDL 1.2.1 path canonicalization rules),
   *  normalize `.` and `..` segments here too.
   */
  private def joinPaths(paths: List[String]): ErrorOr[WomValue] = paths match {
    case Nil => "join_paths requires at least one path".invalidNel
    case first :: rest =>
      rest.find(_.startsWith("/")) match {
        case Some(absolute) =>
          s"join_paths: only the first path may be absolute, but '$absolute' is also absolute".invalidNel
        case None =>
          WomString(rest.foldLeft(first)((joined, next) => joined.stripSuffix("/") + "/" + next)).validNel
      }
  }

  private def readTsv(path: String,
                      hasHeader: Boolean,
                      explicitNames: Option[List[String]],
                      ioFunctionSet: IoFunctionSet,
                      headerElement: ExpressionElement,
                      fieldNamesProvided: Boolean
  ): ErrorOr[WomValue] = {
    val contents: ErrorOr[String] = Try(
      Await.result(
        ioFunctionSet.readFile(path, Option(fileSizeLimitationConfig.readTsvLimit), failOnOverflow = true),
        IoTimeout
      )
    ).toErrorOr.contextualizeErrors(s"""read_tsv("$path")""")

    contents flatMap { text =>
      val lines = text.split("\n").toList.map(_.stripSuffix("\r")).filter(_.nonEmpty)
      val rows = lines.map(_.split("\t", -1).toList)

      (hasHeader, explicitNames) match {
        case (false, None) =>
          headerElement match {
            case PrimitiveLiteralExpressionElement(WomBoolean(false)) if !fieldNamesProvided =>
              WomArray(WomArrayType(WomArrayType(WomStringType)),
                       rows.map(row => WomArray(WomArrayType(WomStringType), row.map(WomString)))
              ).validNel
            case _ =>
              "read_tsv(File, Boolean) requires the Boolean to be true. Use read_tsv(File) to read a file without a header".invalidNel
          }
        case _ =>
          val (fieldNames, dataRows) = (hasHeader, explicitNames) match {
            case (true, Some(names)) => (names, rows.drop(1))
            case (true, None) => (rows.headOption.getOrElse(List.empty), rows.drop(1))
            case (false, names) => (names.getOrElse(List.empty), rows)
          }
          dataRows.zipWithIndex traverse { case (row, index) =>
            if (row.size == fieldNames.size)
              WomObject(fieldNames.zip(row.map(WomString)).toMap).validNel
            else
              s"read_tsv: row ${index + 1} has ${row.size} fields but there are ${fieldNames.size} field names".invalidNel
          } map { objects => WomArray(WomArrayType(WomObjectType), objects) }
      }
    }
  }

  private def sideEffectFile(value: WomValue): Seq[CommandSetupSideEffectFile] = value match {
    case file: WomSingleFile => Seq(CommandSetupSideEffectFile(file))
    case _ => Seq.empty
  }

  /**
   * Writes rows of an `Array[Array[String]]` or an array of structs as TSV, with an optional header row.
   *
   * TODO: Struct fields should be written in the order they are defined, but `WomCompositeType.typeMap` is built from
   *  an unordered `Map` (see `StructEvaluation`), so structs with more than four members may be written in a
   *  different order. Suggestion: keep the struct definition order in `WomCompositeType` (e.g. with a `ListMap`).
   */
  private def writeTsv(array: WomValue,
                       header: Boolean,
                       explicitNames: Option[List[String]],
                       ioFunctionSet: IoFunctionSet
  ): ErrorOr[WomSingleFile] = {
    val rowsAndFieldNames: ErrorOr[(List[List[String]], List[String])] = unwrap(array) match {
      case WomArray(WomArrayType(WomCompositeType(typeMap, _)), structs) =>
        val fieldNames = typeMap.keys.toList
        structs.toList traverse {
          case WomObject(values, _) =>
            fieldNames.map(name => values.get(name).map(unwrap(_).valueString).getOrElse("")).validNel
          case other => s"write_tsv: expected a struct but got ${other.womType.stableName}".invalidNel
        } map (_ -> fieldNames)
      case WomArrayLike(rows) =>
        rows.value.toList traverse {
          case WomArrayLike(row) => row.value.toList.map(unwrap(_).valueString).validNel
          case other => s"write_tsv: expected an array of strings but got ${other.womType.stableName}".invalidNel
        } map (_ -> List.empty)
      case other => s"write_tsv: expected an array but got ${other.womType.stableName}".invalidNel
    }

    rowsAndFieldNames flatMap { case (rows, structFieldNames) =>
      val headerRow: ErrorOr[Option[List[String]]] = (header, explicitNames) match {
        case (false, _) => None.validNel
        case (true, Some(names)) => Option(names).validNel
        case (true, None) if structFieldNames.nonEmpty => Option(structFieldNames).validNel
        case (true, None) =>
          "write_tsv: column names are required to write a header for an Array[Array[String]]".invalidNel
      }
      headerRow flatMap { headerNames =>
        val mismatched = headerNames.toList.flatMap(names => rows.filter(_.size != names.size))
        if (mismatched.nonEmpty)
          s"write_tsv: every row must have ${headerNames.map(_.size).getOrElse(0)} columns to match the header".invalidNel
        else {
          val content = (headerNames.toList ++ rows).map(_.mkString("\t") + "\n").mkString
          import wom.values.HashableString
          Try(Await.result(ioFunctionSet.writeFile(s"write_tsv_${content.md5Sum}.tmp", content), IoTimeout)).toErrorOr
        }
      }
    }
  }

  private def valuesEqual(a: WomValue, b: WomValue): Boolean = (a, b) match {
    case (WomOptionalValue(_, None), WomOptionalValue(_, None)) => true
    case (WomOptionalValue(_, None), _) | (_, WomOptionalValue(_, None)) => false
    case _ =>
      val (left, right) = (unwrap(a), unwrap(b))
      left == right || (left.womType != right.womType && left.valueString == right.valueString && sameKind(left, right))
  }

  // Values of different WOM types that WDL considers equal, e.g. a File and the String of its path
  private def sameKind(left: WomValue, right: WomValue): Boolean = (left, right) match {
    case (_: WomPrimitive, _: WomPrimitive) => true
    case _ => false
  }

  private def containsKey(collection: WomValue, key: WomValue): Boolean = (unwrap(collection), key) match {
    case (WomMap(_, map), k) if map.keys.exists(valuesEqual(_, k)) => true
    case (WomObject(values, _), WomString(name)) => values.contains(name)
    case (nested, WomArrayLike(path)) => containsKeyPath(nested, path.value.map(_.valueString).toList)
    case _ => false
  }

  private def containsKeyPath(collection: WomValue, path: List[String]): Boolean = path match {
    case Nil => true
    case head :: tail =>
      val child: Option[WomValue] = unwrap(collection) match {
        case WomMap(_, map) => map.collectFirst { case (k, v) if k.valueString == head => v }
        case WomObject(values, _) => values.get(head)
        case _ => None
      }
      child match {
        case Some(_) if tail.isEmpty => true
        case Some(WomOptionalValue(_, None)) => false
        case Some(value) => containsKeyPath(value, tail)
        case None => false
      }
  }

  /**
   * Serializes a value for an environment variable. Primitive values use their string form; compound values are
   * serialized as JSON. `None` is exported as an empty string.
   *
   * TODO: The spec refers to "Appendix A: WDL Value Serialization" without saying which serialization applies to
   *  compound values in environment variables. JSON is unambiguous and round-trips, so it is used here. Suggestion:
   *  raise this with the OpenWDL community and align with the other engines' behaviour.
   */
  private def serializeForEnvironment(value: WomValue): String = value match {
    case WomOptionalValue(_, None) => ""
    case WomOptionalValue(_, Some(inner)) => serializeForEnvironment(inner)
    case primitive: WomPrimitive => primitive.valueString
    case compound => ValueEvaluation.valueToJson(compound).compactPrint
  }

  // Single-quote a string for the shell, so that it is never interpreted
  private def shellQuote(value: String): String = "'" + value.replace("'", "'\"'\"'") + "'"

  /**
   * TODO: Computing the size of a `Directory` requires listing its contents, which `IoFunctionSet` does not support
   *  yet. Suggestion: add `IoFunctionSet.directorySize(path)` (or a recursive listing) and implement it in the engine
   *  and backend IO function sets, using the same IO actor commands as directory localization.
   */
  private def filesToSize(value: WomValue): ErrorOr[List[WomFile]] = unwrap(value) match {
    // A String argument is interpreted as the path of a file, as in WDL 1.0 and 1.1
    case s: WomString => List(WomSingleFile(s.value)).validNel
    case other => filesIn(other)
  }

  private def filesIn(value: WomValue): ErrorOr[List[WomFile]] = value match {
    case WomOptionalValue(_, None) => List.empty.validNel
    case WomOptionalValue(_, Some(inner)) => filesIn(inner)
    case d: WomUnlistedDirectory =>
      s"size() of a Directory ('${d.value}') is not supported yet".invalidNel
    case f: WomSingleFile => List(f).validNel
    case WomArrayLike(array) => array.value.toList.flatTraverse(filesIn)
    case WomMap(_, map) => map.toList.flatTraverse { case (k, v) => (filesInKey(k), filesIn(v)).mapN(_ ++ _) }
    case WomPair(left, right) => (filesIn(left), filesIn(right)).mapN(_ ++ _)
    case WomObject(values, _) => values.values.toList.flatTraverse(filesIn)
    case _ => List.empty.validNel
  }

  // Map keys are primitive: only File keys refer to files
  private def filesInKey(key: WomValue): ErrorOr[List[WomFile]] = key match {
    case f: WomSingleFile => List(f).validNel
    case _ => List.empty.validNel
  }
}
