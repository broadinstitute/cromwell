package wdl.types

import scala.collection.JavaConverters._
import wdl.{WdlExpression, WdlSyntaxErrorFormatter}
import wdl.AstTools.EnhancedAstNode
import wom.core.WorkflowSource
import wom.types.{WdlBooleanType, WdlFloatType, WdlIntegerType, WdlType}
import wom.types.WdlType.parser
import wom.values.{WdlBoolean, WdlFloat, WdlInteger, WdlValue}

object WdlFlavoredWomType {
  implicit class FromString(val wdlType: WdlType) extends AnyVal {
    /**
      * Converts WDL source into a WdlValue of this type, if possible.
      *
      * @param workflowSource source code representing the WdlValue
      * @return The WdlValue
      */
    //TODO: return a Try ?
    def fromWorkflowSource(workflowSource: WorkflowSource): WdlValue = {
      wdlType match {
        case WdlFloatType => WdlFloat(workflowSource.toDouble)
        case WdlIntegerType => WdlInteger(workflowSource.toInt)
        case WdlExpressionType => WdlExpression.fromString(workflowSource)
        case WdlBooleanType => WdlBoolean(workflowSource.toBoolean)
        case WdlNamespaceType => ??? // This is what the original code was doing and clearly this is right.
        case _ =>
          val tokens = WdlType.parser.lex(workflowSource, "string")
          val terminalMap = tokens.asScala.toVector.map {(_, workflowSource)}.toMap
          val wdlSyntaxErrorFormatter = WdlSyntaxErrorFormatter(terminalMap)

          /* Parsing as an expression is not sufficient... only a subset of these
           * ASTs are valid as WdlValues and this distinction is done in the
           * .wdlValue() method.
           */
          val ast = WdlType.parser.parse_e(tokens, wdlSyntaxErrorFormatter).toAst

          ast.wdlValue(wdlType, wdlSyntaxErrorFormatter)
      }
    }
  }

  def fromDisplayString(wdlString: String): WdlType = {
    wdlString match {
      case "Expression" => WdlExpressionType
      case _ =>
        val tokens = parser.lex(wdlString, "string")
        val terminalMap = tokens.asScala.toVector.map {(_, wdlString)}.toMap
        val wdlSyntaxErrorFormatter = WdlSyntaxErrorFormatter(terminalMap)

        /* parse_type_e() is the parse function for the $type_e nonterminal in grammar.hgr */
        val ast = parser.parse_type_e(tokens, wdlSyntaxErrorFormatter).toAst

        ast.wdlType(wdlSyntaxErrorFormatter)
    }
  }
}
