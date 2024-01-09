package wdl.draft3.transforms.ast2wdlom

import wdl.draft3.parser.WdlParser.{Ast, AstList, AstNode, Terminal}
import wdl.transforms.base.ast2wdlom.{GenericAst, GenericAstList, GenericAstNode, GenericTerminal}
import scala.jdk.CollectionConverters._

case class Draft3GenericAst(ast: Ast) extends GenericAst {
  override def getAttribute(attr: String): GenericAstNode =
    Option(ast.getAttribute(attr)).map(Draft3GenericAstNode.apply).orNull
  override def getAttributes: Map[String, GenericAstNode] = ast.getAttributes.asScala.toMap collect {
    case (key, value) if value != null => key -> Draft3GenericAstNode(value)
  }
  override def getName: String = ast.getName
}

case class Draft3GenericTerminal(terminal: Terminal) extends GenericTerminal {
  override def getSourceString: String = terminal.getSourceString
  override def getTerminalStr: String = terminal.getTerminalStr
  override def getLine: Int = terminal.getLine
  override def getColumn: Int = terminal.getColumn
}

case class Draft3GenericAstList(astList: AstList) extends GenericAstList {
  override def astNodeList: Seq[GenericAstNode] = astList.asScala.toSeq map Draft3GenericAstNode.apply
}

object Draft3GenericAstNode {
  def apply(astNode: AstNode): GenericAstNode = astNode match {
    case list: AstList => Draft3GenericAstList(list)
    case ast: Ast => Draft3GenericAst(ast)
    case terminal: Terminal => Draft3GenericTerminal(terminal)
    case null =>
      throw new Exception
  }
}
