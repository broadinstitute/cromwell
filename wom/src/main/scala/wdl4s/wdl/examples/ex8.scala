package wdl4s.wdl.examples

import java.nio.file.Paths

import wdl4s.wdl.AstTools
import wdl4s.wdl.AstTools.EnhancedAstNode

object ex8 {
  def main(args: Array[String]): Unit = {
    /* Create syntax tree from contents of file */
    val ast = AstTools.getAst(Paths.get(args(0)))

    /* Second parameter is a descriptor about where the first string came from.
     * Most of the time this would be the URI of where the text was loaded from,
     * but there are no restrictions on what the string can be.
     */
    AstTools.getAst("workflow simple {}", "string")

    /* Print the AST */
    println(ast.toPrettyString)

    /* Traverse the tree to find all Task definitions */
    AstTools.findAsts(ast, "Task") foreach {ast =>
      println(s"Task name: ${ast.getAttribute("name").sourceString}")
    }
  }
}
