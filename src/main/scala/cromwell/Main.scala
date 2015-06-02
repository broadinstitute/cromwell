package cromwell

import java.io.File

import cromwell.binding._
import cromwell.engine.SingleWorkflowRunner
import cromwell.parser.WdlParser.SyntaxError
import cromwell.server.CromwellServer
import spray.json._
import scala.util.{Failure, Success}

object Actions extends Enumeration {
  val parse, validate, run, inputs, server = Value
}

object Main extends App {

  getAction(args.headOption) match {
    case Some(x) if x == Actions.validate => validate(args.tail)
    case Some(x) if x == Actions.inputs => inputs(args.tail)
    case Some(x) if x == Actions.run => run(args.tail)
    case Some(x) if x == Actions.parse => parse(args.tail)
    case Some(x) if x == Actions.server => CromwellServer
    case None => CromwellServer
    case _ => usageAndExit()
  }

  def validate(args: Array[String]): Unit = {
    if (args.length != 1) usageAndExit()
    try {
      val binding = WdlBinding.process(new File(args(0)))
      println("WDL File is valid")
    } catch {
      case e:SyntaxError => println(e)
    }
  }

  def inputs(args: Array[String]): Unit = {
    if (args.length != 1) usageAndExit()
    try {
      import cromwell.binding.types.WdlTypeJsonFormatter._
      val binding = WdlBinding.process(new File(args(0)))
      println(binding.workflows.head.inputs.toJson.prettyPrint)
    } catch {
      case e:SyntaxError => println(e)
    }
  }

  def run(args: Array[String]): Unit = {
    if (args.length != 2) usageAndExit()
    val runner = new SingleWorkflowRunner
    val outputs = runner.run(new File(args(0)), new File(args(1)))
    outputs match {
      case Success(o) =>
        import cromwell.binding.values.WdlValueJsonFormatter._
        println("Workflow Completed.  Outputs are:")
        println(o.toJson.prettyPrint)
      case Failure(f) => println(f.printStackTrace())
    }
  }

  def parse(args: Array[String]): Unit = {
    if (args.length != 1) usageAndExit()
    else println(WdlBinding.getAst(new File(args(0))).toPrettyString)
  }

  def usageAndExit(): Unit = {
    println(
      """
        |java -jar cromwell.jar <action> <parameters>
        |
        |Actions:
        |
        |validate <WDL file>
        |
        |  Performs full validation of the WDL file including syntax
        |  and semantic checking
        |
        |inputs <WDL file>
        |
        |  Print a JSON skeleton file of the inputs needed for this
        |  workflow.  Fill in the values in this JSON document and
        |  pass it in to the 'run' subcommand.
        |
        |run <WDL file> <JSON inputs file>
        |
        |  Given a WDL file and JSON file containing the value of the
        |  workflow inputs, this will run the workflow locally and
        |  print out the outputs in JSON format.
        |
        |parse <WDL file>
        |
        |  Compares a WDL file against the grammar and prints out an
        |  abstract syntax tree if it is valid, and a syntax error
        |  otherwise.  Note that higher-level AST checks are not done
        |  via this sub-command and the 'validate' subcommand should
        |  be used for full validation
        |
        |server
        |
        |  Starts a web server on port 8000.  See the web server
        |  documentation for more details about the API endpoints.
      """.stripMargin)
    System.exit(-1)
  }

  def getAction(firstArg: Option[String]): Option[Actions.Value] = for {
    arg <- firstArg
    a <- Actions.values find { _.toString == arg }
  } yield a
}
