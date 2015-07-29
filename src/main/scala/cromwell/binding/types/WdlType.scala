package cromwell.binding.types

import cromwell.binding.AstTools.EnhancedAstNode
import cromwell.binding.values.WdlValue
import cromwell.binding.{AstTools, WdlSource, WdlSyntaxErrorFormatter}
import cromwell.parser.WdlParser

import scala.collection.JavaConverters._
import scala.util.{Failure, Try}

trait WdlType {

  /**
   * Method to be overridden by implementation classes defining a partial function
   * for the conversion of raw input values to specific implementation class value types.
   * i.e.  `WdlBooleanType` should define a partial function that knows how to
   * construct `WdlBoolean`s for inputs of supported types and contents.  Values for which
   * the partial function is not defined are assumed to not be convertible to the target type.
   */
  protected def coercion: PartialFunction[Any, WdlValue]

  /**
   * Public interface for a `Try`-wrapped conversion of an input of type `Any` to
   * a `WdlValue`.
   */
  def coerceRawValue(any: Any): Try[WdlValue] = {
    if (!coercion.isDefinedAt(any)) {
      Failure(new IllegalArgumentException(s"No coercion defined from '$any' of type '${any.getClass}' to ${getClass.getSimpleName}."))
    } else {
      Try(coercion(any))
    }
  }

  def toWdlString: String

  /**
   * Converts WDL source into a WdlValue of this type, if possible.
   *
   * @param wdlSource source code representing the WdlValue
   * @return The WdlValue
   */
  def fromWdlString(wdlSource: WdlSource): WdlValue = {
    val tokens = WdlType.parser.lex(wdlSource, "string")
    val terminalMap = tokens.asScala.toVector.map {(_, wdlSource)}.toMap
    val wdlSyntaxErrorFormatter = new WdlSyntaxErrorFormatter(terminalMap)

    /* Parsing as an expression is not sufficient... only a subset of these
     * ASTs are valid as WdlValues and this distinction is done in the
     * .wdlValue() method.
     */
    val ast = WdlType.parser.parse_e(tokens, wdlSyntaxErrorFormatter).toAst

    ast.wdlValue(this, wdlSyntaxErrorFormatter)
  }
}

object WdlType {
  val parser = new WdlParser()

  private lazy val wdlTypes = Seq(
    WdlBooleanType, WdlExpressionType, WdlFileType, WdlFloatType,
    WdlIntegerType, WdlNamespaceType, WdlObjectType, WdlStringType
  )

  def fromWdlString(wdlString: String): WdlType = {
    wdlString match {
      case "Expression" => WdlExpressionType
      case _ =>
        val tokens = parser.lex(wdlString, "string")
        val terminalMap = tokens.asScala.toVector.map {(_, wdlString)}.toMap
        val wdlSyntaxErrorFormatter = new WdlSyntaxErrorFormatter(terminalMap)

        /* parse_type_e() is the parse function for the $type_e nonterminal in grammar.hgr */
        val ast = parser.parse_type_e(tokens, wdlSyntaxErrorFormatter).toAst

        ast.wdlType(wdlSyntaxErrorFormatter)
    }
  }
}
