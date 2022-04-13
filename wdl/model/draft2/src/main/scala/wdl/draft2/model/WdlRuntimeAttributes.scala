package wdl.draft2.model

import common.collections.EnhancedCollections._
import wdl.draft2.model.AstTools.{AstNodeName, EnhancedAstNode}
import wdl.draft2.parser.WdlParser.{Ast, AstList}
import wom.RuntimeAttributes

import scala.jdk.CollectionConverters._

case class WdlRuntimeAttributes(attrs: Map[String, WdlExpression]) {
  def toWomRuntimeAttributes(task: WdlTask) = RuntimeAttributes(attrs.safeMapValues(WdlWomExpression(_, task)))
}

object WdlRuntimeAttributes {
  def apply(ast: Ast): WdlRuntimeAttributes = {
    val asts = ast.findAsts(AstNodeName.Runtime)
    if (asts.size > 1) throw new UnsupportedOperationException("Only one runtime block may be defined per task")
    val kvPairAsts = asts.headOption.map(_.getAttribute("map").asInstanceOf[AstList].asScala.toVector.map(_.asInstanceOf[Ast]))
    val runtimeAttributeMap: Map[String, WdlExpression] = kvPairAsts match {
      case Some(vector) => vector.map(ast => processRuntimeAttribute(ast)).toMap
      case None => Map.empty[String, WdlExpression]
    }

    WdlRuntimeAttributes(runtimeAttributeMap)
  }

  private def processRuntimeAttribute(ast: Ast): (String, WdlExpression) = {
    val key = ast.getAttribute("key").sourceString
    val expression = WdlExpression(ast.getAttribute("value"))
    key -> expression
  }
}
