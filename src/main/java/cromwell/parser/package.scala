package cromwell

/**
 * ==Workflow Description Language Parser==
 *
 * This package contains the parser for WDL, which is auto-generated from the ''grammar.hgr'' file in ''src/main/java/cromwell/resources''.
 *
 * The main interface to this package are the following functions:
 *
 * - lex() which takes WDL source and tokenizes it
 * - parse() which takes the output of lex(), a sequence of tokens, and returns a parse tree
 *
 * Here's a quick example of using this package to obtain an Abstract Syntax Tree of a WDL file:
 *
 * {{{
 * val wdlFile = scala.io.Source.fromFile("/path/to/file.wdl")
 * val wdlContents = try wdlFile.mkString finally wdlFile.close()
 * val parser = new WdlParser()
 * val tokens = new WdlParser.TokenStream(parser.lex(wdlContents, args(1)))
 * val ast = parser.parse(tokens).toAst()
 * println(ast.toPrettyString())
 * }}}
 *
 */

package object parser {

}
