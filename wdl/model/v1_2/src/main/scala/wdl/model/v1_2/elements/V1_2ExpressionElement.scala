package wdl.model.v1_2.elements

import wdl.model.draft3.elements.ExpressionElement
import wdl.shared.model.expression.ExtensionExpressionElement

/**
 * Expression elements introduced in WDL 1.2 (https://github.com/openwdl/wdl/releases/tag/v1.2.1).
 *
 * The element hierarchy in `wdl.model.draft3.elements` is shared with WDL 1.0 and 1.1 and is intentionally left
 * unchanged, so elements that only exist in WDL 1.2 are defined here.
 */
object V1_2ExpressionElement {

  /** A WDL 1.2 operator or function call. `params` are its argument expressions, in order. */
  sealed trait V1_2FunctionCallElement extends ExpressionElement with ExtensionExpressionElement[ExpressionElement] {

    /** The WDL function name, used to write the element as WDL source. */
    def functionName: String

    override def toWdlV1(paramToWdl: ExpressionElement => String): String =
      s"$functionName(${params.map(paramToWdl).mkString(", ")})"
  }

  // Operators

  /** `x ** y` */
  final case class Exponentiation(left: ExpressionElement, right: ExpressionElement) extends V1_2FunctionCallElement {
    override def functionName: String = "**"
    override def params: Seq[ExpressionElement] = Seq(left, right)
    override def withParams(newParams: Seq[ExpressionElement]): ExpressionElement =
      Exponentiation(newParams(0), newParams(1))
    override def toWdlV1(paramToWdl: ExpressionElement => String): String =
      s"(${paramToWdl(left)} ** ${paramToWdl(right)})"
  }

  // String functions

  /** `String? find(String input, String pattern)` */
  final case class Find(input: ExpressionElement, pattern: ExpressionElement) extends V1_2FunctionCallElement {
    override def functionName: String = "find"
    override def params: Seq[ExpressionElement] = Seq(input, pattern)
    override def withParams(newParams: Seq[ExpressionElement]): ExpressionElement = Find(newParams(0), newParams(1))
  }

  /** `Boolean matches(String input, String pattern)` */
  final case class Matches(input: ExpressionElement, pattern: ExpressionElement) extends V1_2FunctionCallElement {
    override def functionName: String = "matches"
    override def params: Seq[ExpressionElement] = Seq(input, pattern)
    override def withParams(newParams: Seq[ExpressionElement]): ExpressionElement = Matches(newParams(0), newParams(1))
  }

  // File functions

  /**
   * `String join_paths(Directory, String)`, `String join_paths(Directory, Array[String]+)` or
   * `String join_paths(Array[String]+)`.
   */
  final case class JoinPaths(first: ExpressionElement, rest: Option[ExpressionElement])
      extends V1_2FunctionCallElement {
    override def functionName: String = "join_paths"
    override def params: Seq[ExpressionElement] = first +: rest.toSeq
    override def withParams(newParams: Seq[ExpressionElement]): ExpressionElement =
      JoinPaths(newParams.head, newParams.lift(1))
  }

  /** `Array[Object] read_tsv(File, Boolean[, Array[String]])` */
  final case class ReadTsvWithHeader(file: ExpressionElement,
                                     header: ExpressionElement,
                                     fieldNames: Option[ExpressionElement]
  ) extends V1_2FunctionCallElement {
    override def functionName: String = "read_tsv"
    override def params: Seq[ExpressionElement] = Seq(file, header) ++ fieldNames.toSeq
    override def withParams(newParams: Seq[ExpressionElement]): ExpressionElement =
      ReadTsvWithHeader(newParams(0), newParams(1), newParams.lift(2))
  }

  /** `File write_tsv(Array[Array[String]]|Array[Struct], Boolean[, Array[String]])` */
  final case class WriteTsvWithHeader(array: ExpressionElement,
                                      header: ExpressionElement,
                                      columnNames: Option[ExpressionElement]
  ) extends V1_2FunctionCallElement {
    override def functionName: String = "write_tsv"
    override def params: Seq[ExpressionElement] = Seq(array, header) ++ columnNames.toSeq
    override def withParams(newParams: Seq[ExpressionElement]): ExpressionElement =
      WriteTsvWithHeader(newParams(0), newParams(1), newParams.lift(2))
  }

  // Array functions

  /** `Boolean contains(Array[P], P)` */
  final case class Contains(array: ExpressionElement, value: ExpressionElement) extends V1_2FunctionCallElement {
    override def functionName: String = "contains"
    override def params: Seq[ExpressionElement] = Seq(array, value)
    override def withParams(newParams: Seq[ExpressionElement]): ExpressionElement = Contains(newParams(0), newParams(1))
  }

  /** `Array[Array[X]] chunk(Array[X], Int)` */
  final case class Chunk(array: ExpressionElement, size: ExpressionElement) extends V1_2FunctionCallElement {
    override def functionName: String = "chunk"
    override def params: Seq[ExpressionElement] = Seq(array, size)
    override def withParams(newParams: Seq[ExpressionElement]): ExpressionElement = Chunk(newParams(0), newParams(1))
  }

  /** `X select_first(Array[X?], X)` */
  final case class SelectFirstWithDefault(array: ExpressionElement, default: ExpressionElement)
      extends V1_2FunctionCallElement {
    override def functionName: String = "select_first"
    override def params: Seq[ExpressionElement] = Seq(array, default)
    override def withParams(newParams: Seq[ExpressionElement]): ExpressionElement =
      SelectFirstWithDefault(newParams(0), newParams(1))
  }

  // Map functions

  /** `Boolean contains_key(Map[P, Y]|Object|Struct, P|Array[String])` */
  final case class ContainsKey(collection: ExpressionElement, key: ExpressionElement) extends V1_2FunctionCallElement {
    override def functionName: String = "contains_key"
    override def params: Seq[ExpressionElement] = Seq(collection, key)
    override def withParams(newParams: Seq[ExpressionElement]): ExpressionElement =
      ContainsKey(newParams(0), newParams(1))
  }

  /** `Array[Y] values(Map[P, Y])` */
  final case class Values(map: ExpressionElement) extends V1_2FunctionCallElement {
    override def functionName: String = "values"
    override def params: Seq[ExpressionElement] = Seq(map)
    override def withParams(newParams: Seq[ExpressionElement]): ExpressionElement = Values(newParams(0))
  }

  // Engine-internal elements (not callable from WDL source)

  /**
   * The value of a declaration annotated with the WDL 1.2 `env` modifier, serialized and quoted so that it can be
   * safely assigned to an environment variable in the task's command script.
   */
  final case class EnvironmentVariableValue(value: ExpressionElement) extends V1_2FunctionCallElement {
    // Only used to describe the element, e.g. in the command template. It cannot be called from WDL source.
    override def functionName: String = "__environment_variable_value"
    override def params: Seq[ExpressionElement] = Seq(value)
    override def withParams(newParams: Seq[ExpressionElement]): ExpressionElement =
      EnvironmentVariableValue(newParams(0))
  }
}
