package cromwell

import cromwell.parser.WdlParser

object Main {
    def main(args: Array[String]) {
        val actions = List("parse", "run")
        if (args.length < 2 || !actions.contains(args(0))) {
            println("Usage: cromwell.jar parse <wdl file>")
            println("Usage: cromwell.jar run <wdl file>")
            System.exit(-1)
        }
        if (args(0).equals("parse")) {
            val wdlFile = scala.io.Source.fromFile(args(1))
            val wdlContents = try wdlFile.mkString finally wdlFile.close()
            val parser = new WdlParser()
            val tokens = new WdlParser.TokenStream(parser.lex(wdlContents, args(1)))
            val ast = parser.parse(tokens).toAst()
            println(ast.toPrettyString())
        }
        if (args(0).equals("run")) {

        }
    }
}
