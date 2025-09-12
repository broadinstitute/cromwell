package wdl.draft2.model

import common.collections.EnhancedCollections._
import wdl.draft2.model.AstTools.{AstNodeName, EnhancedAstNode}
import wdl.draft2.parser.WdlParser.{Ast, AstList}
import wom.{RuntimeAttributes, RuntimeAttributesKeys}

import scala.jdk.CollectionConverters._

case class WdlRuntimeAttributes(attrs: Map[String, WdlExpression]) {
  def toWomRuntimeAttributes(task: WdlTask) =
    // In future WDL versions, `container` has superceded `docker` as the runtime attribute for specifying a
    // container image. For pre-1.1 WDLs, we need to remove `container` if it exists so that it doesn't interfere
    // with `docker`.
    RuntimeAttributes(
      attrs.filterNot(m => m._1 == RuntimeAttributesKeys.ContainerKey).safeMapValues(WdlWomExpression(_, task))
    )
}

object WdlRuntimeAttributes {
  def apply(ast: Ast): WdlRuntimeAttributes = {
    val asts = ast.findAsts(AstNodeName.Runtime)
    if (asts.size > 1) throw new UnsupportedOperationException("Only one runtime block may be defined per task")
    val kvPairAsts =
      asts.headOption.map(_.getAttribute("map").asInstanceOf[AstList].asScala.toVector.map(_.asInstanceOf[Ast]))
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
