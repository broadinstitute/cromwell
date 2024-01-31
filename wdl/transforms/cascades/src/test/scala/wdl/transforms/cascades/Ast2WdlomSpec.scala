package wdl.transforms.cascades

import java.util

import common.Checked
import common.assertion.ErrorOrAssertions._
import common.transforms.CheckedAtoB
import cats.instances.either._
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.cascades.parser.WdlParser
import wdl.cascades.parser.WdlParser.{ParseTree, SyntaxErrorFormatter}
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.elements._
import wdl.transforms.base.ast2wdlom.{GenericAstNode, _}
import wdl.transforms.cascades.Ast2WdlomSpec._
import wdl.transforms.cascades.ast2wdlom._
import wdl.transforms.cascades.parsing.WdlCascadesSyntaxErrorFormatter
import wom.callable.MetaValueElement.MetaValueElementInteger
import wom.types.WomIntegerType
import wom.values.WomInteger

import scala.jdk.CollectionConverters._

object Ast2WdlomSpec {
  val parser = new WdlParser()

  def fromString[A](expression: String,
                    parseFunction: (util.List[WdlParser.Terminal], SyntaxErrorFormatter) => ParseTree
  )(implicit converter: CheckedAtoB[GenericAstNode, A]): Checked[A] = {

    // Add the "version development" to force the lexer into "main" mode.
    val versionedExpression = "version development\n" + expression
    // That "version development" means we'll have 2 unwanted tokens at the start of the list, so drop 'em:
    val tokens = parser.lex(versionedExpression, "string").asScala.drop(2).asJava
    val terminalMap = (tokens.asScala.toVector map { (_, versionedExpression) }).toMap
    val parseTree = parseFunction(tokens, WdlCascadesSyntaxErrorFormatter(terminalMap))
    (wrapAstNode andThen converter).run(parseTree.toAst)
  }
}

class Ast2WdlomSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  it should "parse a simple expression" in {
    val str = "3 + 3"
    val expr = fromString[ExpressionElement](str, parser.parse_e)
    expr shouldBeValid Add(PrimitiveLiteralExpressionElement(WomInteger(3)),
                           PrimitiveLiteralExpressionElement(WomInteger(3))
    )
  }

  it should "parse a map expression" in {
    val str = "{3: 3}"
    val expr = fromString[ExpressionElement](str, parser.parse_e)
    expr shouldBeValid MapLiteral(
      Map(PrimitiveLiteralExpressionElement(WomInteger(3)) -> PrimitiveLiteralExpressionElement(WomInteger(3)))
    )
  }

  it should "parse a simple meta section" in {
    val str = "meta { five: 5 }"
    val metaKvPair = fromString[MetaSectionElement](str, parser.parse_meta)
    metaKvPair shouldBeValid MetaSectionElement(Map("five" -> MetaValueElementInteger(5)))
  }

  it should "parse a struct element" in {
    val str = "struct Foo { Int five\n Int six }"
    val struct = fromString[StructElement](str, parser.parse_struct)(astNodeToAst andThen astToStructElement)
    struct shouldBeValid StructElement("Foo",
                                       Seq(StructEntryElement("five", PrimitiveTypeElement(WomIntegerType)),
                                           StructEntryElement("six", PrimitiveTypeElement(WomIntegerType))
                                       )
    )
  }

  it should "parse the new as_map function" in {
    val str = "as_map(some_pairs)"
    val expr = fromString[ExpressionElement](str, parser.parse_e)
    expr shouldBeValid (AsMap(IdentifierLookup("some_pairs")))
  }

  it should "parse the new as_pairs function" in {
    val str = "as_pairs(some_map)"
    val expr = fromString[ExpressionElement](str, parser.parse_e)
    expr shouldBeValid (AsPairs(IdentifierLookup("some_map")))
  }

  it should "parse the new collect_by_key function" in {
    val str = "collect_by_key(some_map)"
    val expr = fromString[ExpressionElement](str, parser.parse_e)
    expr shouldBeValid (CollectByKey(IdentifierLookup("some_map")))
  }

  it should "parse the new None keyword" in {
    val str = "None"
    val expr = fromString[ExpressionElement](str, parser.parse_e)
    expr shouldBeValid NoneLiteralElement
  }

  it should "parse the new suffix function" in {
    val str = "suffix(some_str, some_arr)"
    val expr = fromString[ExpressionElement](str, parser.parse_e)
    expr shouldBeValid (Suffix(IdentifierLookup("some_str"), IdentifierLookup("some_arr")))
  }
}
