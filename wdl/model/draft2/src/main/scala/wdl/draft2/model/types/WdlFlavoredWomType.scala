package wdl.draft2.model.types

import wdl.draft2.model.AstTools.EnhancedAstNode
import wdl.draft2.model.{WdlExpression, WdlSyntaxErrorFormatter}
import wdl.draft2.parser.WdlParser
import wom.core.WorkflowSource
import wom.types.{WomBooleanType, WomFloatType, WomIntegerType, WomType}
import wom.values.{WomBoolean, WomFloat, WomInteger, WomValue}

import scala.collection.JavaConverters._

object WdlFlavoredWomType {
  private val parser = new WdlParser()

  implicit class FromString(val womType: WomType) extends AnyVal {
    /**
      * Converts WDL source into a WomValue of this type, if possible.
      *
      * @param workflowSource source code representing the WomValue
      * @return The WomValue
      */
    //TODO: return a Try ?
    def fromWorkflowSource(workflowSource: WorkflowSource): WomValue = {
      womType match {
        case WomFloatType => WomFloat(workflowSource.toDouble)
        case WomIntegerType => WomInteger(workflowSource.toInt)
        case WdlExpressionType => WdlExpression.fromString(workflowSource)
        case WomBooleanType => WomBoolean(workflowSource.toBoolean)
        case WdlNamespaceType => throw new UnsupportedOperationException // This is what the original code was doing and clearly this is right.
        case _ =>
          val tokens = parser.lex(workflowSource, "string")
          val terminalMap = tokens.asScala.toVector.map {(_, workflowSource)}.toMap
          val wdlSyntaxErrorFormatter = WdlSyntaxErrorFormatter(terminalMap)

          /* Parsing as an expression is not sufficient... only a subset of these
           * ASTs are valid as WdlValues and this distinction is done in the
           * .womValue() method.
           */
          val ast = parser.parse_e(tokens, wdlSyntaxErrorFormatter).toAst

          ast.womValue(womType, wdlSyntaxErrorFormatter)
      }
    }
  }

  def fromDisplayString(wdlString: String): WomType = {
    wdlString match {
      case "Expression" => WdlExpressionType
      case _ =>
        val tokens = parser.lex(wdlString, "string")
        val terminalMap = tokens.asScala.toVector.map {(_, wdlString)}.toMap
        val wdlSyntaxErrorFormatter = WdlSyntaxErrorFormatter(terminalMap)

        /* parse_type_e() is the parse function for the $type_e nonterminal in grammar.hgr */
        val ast = parser.parse_type_e(tokens, wdlSyntaxErrorFormatter).toAst

        ast.womType(wdlSyntaxErrorFormatter)
    }
  }
}
