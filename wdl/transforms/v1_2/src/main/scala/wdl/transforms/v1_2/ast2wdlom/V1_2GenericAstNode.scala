package wdl.transforms.v1_2.ast2wdlom

import wdl.v1_2.parser.WdlParser.{Ast, AstList, AstNode, Terminal}
import wdl.transforms.base.ast2wdlom.{GenericAst, GenericAstList, GenericAstNode, GenericTerminal}
import scala.jdk.CollectionConverters._

case class V1_2GenericAst(ast: Ast) extends GenericAst {
  override def getAttribute(attr: String): GenericAstNode =
    Option(ast.getAttribute(attr)).map(V1_2GenericAstNode.apply).orNull
  override def getAttributes: Map[String, GenericAstNode] = ast.getAttributes.asScala.toMap collect {
    case (key, value) if value != null => key -> V1_2GenericAstNode(value)
  }
  override def getName: String = ast.getName
}

case class V1_2GenericTerminal(terminal: Terminal) extends GenericTerminal {
  override def getSourceString: String = terminal.getSourceString
  // `task` is a keyword, but in an expression it refers to the hidden `task` variable (WDL 1.2)
  override def getTerminalStr: String = terminal.getTerminalStr match {
    case "task" => "identifier"
    case other => other
  }
  override def getLine: Int = terminal.getLine
  override def getColumn: Int = terminal.getColumn

  override def toString: String = s"V1_2GenericTerminal($getTerminalStr: $getSourceString)"
}

case class V1_2GenericAstList(astList: AstList) extends GenericAstList {
  override def astNodeList: Seq[GenericAstNode] = astList.asScala.toSeq map V1_2GenericAstNode.apply
}

/**
 * Nodes that do not exist in the Hermes parse tree. They let new WDL 1.2 syntax be presented to the shared
 * (WDL 1.0+) AST converters in a shape those converters already understand, without changing the shared code.
 */
final case class SyntheticGenericAst(name: String, attributes: Map[String, GenericAstNode]) extends GenericAst {
  override def getAttribute(attr: String): GenericAstNode = attributes.get(attr).orNull
  override def getAttributes: Map[String, GenericAstNode] = attributes
  override def getName: String = name
}

final case class SyntheticGenericAstList(astNodeList: Seq[GenericAstNode]) extends GenericAstList

final case class SyntheticGenericTerminal(getTerminalStr: String, getSourceString: String, getLine: Int, getColumn: Int)
    extends GenericTerminal

object V1_2GenericAstNode {

  /**
   * Name of the synthetic engine function used to represent the `**` operator. It is not a valid WDL identifier, so it
   * can never clash with a function called from WDL source.
   */
  val ExponentiationFunctionName = "**"

  def apply(astNode: AstNode): GenericAstNode = astNode match {
    case list: AstList => V1_2GenericAstList(list)
    case ast: Ast if ast.getName == "Exponentiation" => exponentiationAsFunctionCall(ast)
    case ast: Ast if ast.getName == "MultilineStringLiteral" => MultilineStrings.toStringLiteral(ast)
    case ast: Ast => V1_2GenericAst(ast)
    case terminal: Terminal => V1_2GenericTerminal(terminal)
    case null =>
      throw new Exception()
  }

  // `x ** y` is presented to the shared expression converter as the function call `**(x, y)`
  private def exponentiationAsFunctionCall(ast: Ast): GenericAst = {
    val lhs = V1_2GenericAstNode(ast.getAttribute("lhs"))
    val rhs = V1_2GenericAstNode(ast.getAttribute("rhs"))
    val position = lhs.firstTerminal
    val name = SyntheticGenericTerminal("identifier",
                                        ExponentiationFunctionName,
                                        position.map(_.getLine).getOrElse(0),
                                        position.map(_.getColumn).getOrElse(0)
    )
    SyntheticGenericAst("FunctionCall", Map("name" -> name, "params" -> SyntheticGenericAstList(Seq(lhs, rhs))))
  }
}
