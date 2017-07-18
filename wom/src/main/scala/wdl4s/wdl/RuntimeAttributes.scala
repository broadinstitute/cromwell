package wdl4s.wdl

import wdl4s.parser.WdlParser.{Ast, AstList}
import wdl4s.wdl.AstTools.{AstNodeName, EnhancedAstNode}

import scala.collection.JavaConverters._

case class RuntimeAttributes(attrs: Map[String, WdlExpression])

object RuntimeAttributes {
  def apply(ast: Ast): RuntimeAttributes = {
    val asts = ast.findAsts(AstNodeName.Runtime)
    if (asts.size > 1) throw new UnsupportedOperationException("Only one runtime block may be defined per task")
    val kvPairAsts = asts.headOption.map(_.getAttribute("map").asInstanceOf[AstList].asScala.toVector.map(_.asInstanceOf[Ast]))
    val runtimeAttributeMap: Map[String, WdlExpression] = kvPairAsts match {
      case Some(vector) => vector.map(ast => processRuntimeAttribute(ast)).toMap
      case None => Map.empty[String, WdlExpression]
    }

    RuntimeAttributes(runtimeAttributeMap)
  }

  private def processRuntimeAttribute(ast: Ast): (String, WdlExpression) = {
    val key = ast.getAttribute("key").sourceString
    val expression = WdlExpression(ast.getAttribute("value"))
    key -> expression
  }
}
