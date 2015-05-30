package cromwell

import java.io.File

import cromwell.binding._
import cromwell.parser.WdlParser.SyntaxError
import cromwell.server.CromwellServer


object Actions extends Enumeration {
  val parse, validate, run, server = Value
}


object Main extends App {

  getAction(args.headOption) match {
    case Some(x) if x == Actions.parse => parseIt(args.tail)
    case Some(x) if x == Actions.validate => validateIt(args.tail)
    case Some(x) if x == Actions.run => runIt(args.tail)
    case Some(x) if x == Actions.server => CromwellServer // Mention it so it gets instantiated
    case None => CromwellServer
    case _ => usageAndExit()
  }

  def parseIt(args: Array[String]): Unit = {
    if (args.isEmpty) usageAndExit()
    else println(WdlBinding.getAst(new File(args(0))).toPrettyString)
  }

  def validateIt(args: Array[String]): Unit = {
    try {
      val binding = WdlBinding.process(new File(args(0)))
      println("WDL File is valid")
    } catch {
      case e:SyntaxError => println(e)
    }
    if (args.isEmpty) usageAndExit()
  }

  def runIt(args: Array[String]): Unit = {
    throw new UnsupportedOperationException("Not yet implemented.")
  }

  def usageAndExit(): Unit = {
    println("Usage: cromwell.jar parse <wdl file>")
    println("Usage: cromwell.jar validate <wdl file>")
    println("Usage: cromwell.jar run <wdl file> <JSON inputs>")
    println("Usage: cromwell.jar server")
    System.exit(-1)
  }

  def getAction(firstArg: Option[String]): Option[Actions.Value] = for {
    arg <- firstArg
    a <- Actions.values find { _.toString == arg }
  } yield a
}
