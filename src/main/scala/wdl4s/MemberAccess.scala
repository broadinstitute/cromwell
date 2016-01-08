package wdl4s

import wdl4s.parser.WdlParser.{Terminal, Ast}
import wdl4s.AstTools.EnhancedAstNode

/**
 * Represents a member-access construct holding the left and right hand sides, e.g.
 * "foo.blah" would deconstruct into "foo" and "blah".
 *
 * The right-hand side of a member-access AST should always be interpreted as a String
 * Sometimes, the left-hand side is itself a MemberAccess AST, like in the expression
 * for `call t1` below.  In that example, "ns.ns2.task_name" would be the left hand side
 * In the `call t2` example, "alias" would be the left hand side
 *
 * import "test.wdl" as ns
 * workflow w {
 *  call ns.ns2.task_name
 *  call t1 {
 *    input: x=ns.ns2.task_name.output
 *  }
 *
 *  call ns.ns2.task_name as alias
 *  call t2 {
 *    input: y=alias.output
 *  }
 *}
 */
case class MemberAccess(lhs: String, rhs: String)

object MemberAccess {
  def apply(ast: Ast): MemberAccess = {
    val rhs = ast.getAttribute("rhs").sourceString

    val lhs = ast.getAttribute("lhs") match {
      case a: Ast => WdlExpression.toString(a)
      case terminal: Terminal => terminal.sourceString
    }

    MemberAccess(lhs, rhs)
  }
}