package cromwell.binding.types

import cromwell.binding.{WdlExpressionException, WdlFunctions, WdlExpression}
import cromwell.binding.values.{WdlInteger, WdlArray}
import cromwell.parser.WdlParser.{AstList, Ast}
import spray.json.JsArray
import scala.collection.JavaConverters._
import scala.util.{Success, Failure}

case class WdlArrayType(memberType: WdlType) extends WdlType {
  val toWdlString: String = s"Array[${memberType.toWdlString}]"

  private def coerceIterable(s: Seq[Any]): WdlArray = {
    val coerced = s.map {memberType.coerceRawValue(_).get}
    WdlArray(WdlArrayType(coerced.head.wdlType), coerced)
  }

  override protected def coercion = {
    case s: Seq[Any] if s.nonEmpty => coerceIterable(s)
    case js: JsArray if js.elements.nonEmpty => coerceIterable(js.elements)
  }

  override def fromRawString(str: String) = {
    val expr: WdlExpression = WdlExpression.fromString(str)
    if (!expr.ast.isInstanceOf[Ast]) throw new WdlExpressionException(s"Could not convert $str to a ${this.toWdlString}")
    if (expr.ast.asInstanceOf[Ast].getName != "ArrayLiteral") throw new WdlExpressionException(s"Invalid AST to convert to ${this.toWdlString}:\n\n${expr.ast.asInstanceOf[Ast].toPrettyString}")
    val ast = expr.ast.asInstanceOf[Ast]
    class NoFunctions extends WdlFunctions {
      def getFunction(name: String): WdlFunction = throw new WdlExpressionException("TODO: array definitions currently do not support function calls")
    }
    val elements = ast.getAttribute("values").asInstanceOf[AstList].asScala.toVector.map {expr =>
      WdlExpression(expr).evaluate((s:String) => throw new WdlExpressionException("TODO: array definitions currently do not support variable lookups"), new NoFunctions())
    }
    WdlArray(this, elements.map {
      case Success(v) => v
      case Failure(f) => throw new WdlExpressionException(s"Failed to evaluate expression: $f")
    })
  }
}
