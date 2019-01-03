package wdl.transforms.biscayne.ast2wdlom

import wdl.biscayne.parser.WdlParser.{Ast, AstList, AstNode, Terminal}
import wdl.transforms.base.ast2wdlom.{GenericAst, GenericAstList, GenericAstNode, GenericTerminal}
import scala.collection.JavaConverters._

case class BiscayneGenericAst(ast: Ast) extends GenericAst {
  override def getAttribute(attr: String): GenericAstNode = Option(ast.getAttribute(attr)).map(BiscayneGenericAstNode.apply).orNull
  override def getAttributes: Map[String, GenericAstNode] = ast.getAttributes.asScala.toMap collect {
    case (key, value) if value != null => key -> BiscayneGenericAstNode(value)
  }
  override def getName: String = ast.getName
}

case class BiscayneGenericTerminal(terminal: Terminal) extends GenericTerminal {
  override def getSourceString: String = terminal.getSourceString
  override def getTerminalStr: String = terminal.getTerminalStr
  override def getLine: Int = terminal.getLine
  override def getColumn: Int = terminal.getColumn

  override def toString: String = s"BiscayneGenericTerminal($getTerminalStr: $getSourceString)"
}

case class BiscayneGenericAstList(astList: AstList) extends GenericAstList {
  override def astNodeList: Seq[GenericAstNode] = astList.asScala.toSeq map BiscayneGenericAstNode.apply
}

object BiscayneGenericAstNode {
  def apply(astNode: AstNode): GenericAstNode = astNode match {
    case list: AstList => BiscayneGenericAstList(list)
    case ast: Ast => BiscayneGenericAst(ast)
    case terminal: Terminal => BiscayneGenericTerminal(terminal)
    case null =>
      throw new Exception()
  }
}
