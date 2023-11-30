package wdl.transforms.cascades.ast2wdlom

import wdl.cascades.parser.WdlParser.{Ast, AstList, AstNode, Terminal}
import wdl.transforms.base.ast2wdlom.{GenericAst, GenericAstList, GenericAstNode, GenericTerminal}
import scala.jdk.CollectionConverters._

case class cascadesGenericAst(ast: Ast) extends GenericAst {
  override def getAttribute(attr: String): GenericAstNode =
    Option(ast.getAttribute(attr)).map(cascadesGenericAstNode.apply).orNull
  override def getAttributes: Map[String, GenericAstNode] = ast.getAttributes.asScala.toMap collect {
    case (key, value) if value != null => key -> cascadesGenericAstNode(value)
  }
  override def getName: String = ast.getName
}

case class cascadesGenericTerminal(terminal: Terminal) extends GenericTerminal {
  override def getSourceString: String = terminal.getSourceString
  override def getTerminalStr: String = terminal.getTerminalStr
  override def getLine: Int = terminal.getLine
  override def getColumn: Int = terminal.getColumn

  override def toString: String = s"cascadesGenericTerminal($getTerminalStr: $getSourceString)"
}

case class cascadesGenericAstList(astList: AstList) extends GenericAstList {
  override def astNodeList: Seq[GenericAstNode] = astList.asScala.toSeq map cascadesGenericAstNode.apply
}

object cascadesGenericAstNode {
  def apply(astNode: AstNode): GenericAstNode = astNode match {
    case list: AstList => cascadesGenericAstList(list)
    case ast: Ast => cascadesGenericAst(ast)
    case terminal: Terminal => cascadesGenericTerminal(terminal)
    case null =>
      throw new Exception()
  }
}
