package wdl.draft3.transforms.ast2wdlom

import java.util

import cats.instances.either._
import common.Checked
import common.assertion.ErrorOrAssertions._
import common.transforms.CheckedAtoB
import org.scalatest.{FlatSpec, Matchers}
import wdl.draft3.parser.WdlParser
import wdl.draft3.parser.WdlParser.{ParseTree, SyntaxErrorFormatter}
import wdl.draft3.transforms.parsing.WdlDraft3SyntaxErrorFormatter
import wdl.model.draft3.elements.ExpressionElement.IdentifierLookup
import wdl.model.draft3.elements._
import wdl.transforms.base.ast2wdlom.GenericAstNode

import scala.collection.JavaConverters._

class Ast2WdlomSpec extends FlatSpec with Matchers {

  val parser = new WdlParser()

  def fromString[A](expression: String,
                    parseFunction: (util.List[WdlParser.Terminal], SyntaxErrorFormatter) => ParseTree)
                   (implicit converter: CheckedAtoB[GenericAstNode, A]): Checked[A] = {
    // Add the "version 1.0" to force the lexer into "main" mode.
    val versionedExpression = "version 1.0\n" + expression
    // That "version 1.0" means we'll have 2 unwanted tokens at the start of the list, so drop 'em:
    val tokens = parser.lex(versionedExpression, "string").asScala.drop(2).asJava
    val terminalMap = (tokens.asScala.toVector map {(_, expression)}).toMap
    val parseTree = parseFunction(tokens, WdlDraft3SyntaxErrorFormatter(terminalMap))
    (wrapAstNode andThen converter).run(parseTree.toAst)
  }

  it should "not parse the new as_map function" in {
    val str = "as_map(some_pairs)"
    val expr = fromString[ExpressionElement](str, parser.parse_e)
    expr shouldBeInvalid "Failed to parse expression (reason 1 of 1): Unknown engine function: 'as_map'"
  }

  it should "not parse the new as_pairs function" in {
    val str = "as_pairs(some_map)"
    val expr = fromString[ExpressionElement](str, parser.parse_e)
    expr shouldBeInvalid "Failed to parse expression (reason 1 of 1): Unknown engine function: 'as_pairs'"
  }

  it should "not parse the new collect_by_key function" in {
    val str = "collect_by_key(some_map)"
    val expr = fromString[ExpressionElement](str, parser.parse_e)
    expr shouldBeInvalid "Failed to parse expression (reason 1 of 1): Unknown engine function: 'collect_by_key'"
  }

  it should "parse the (biscayne) None keyword as a plain old identifier" in {
    val str = "None"
    val expr = fromString[ExpressionElement](str, parser.parse_e)
    expr shouldBeValid(IdentifierLookup("None"))
  }
}
