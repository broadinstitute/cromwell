package cromwell

import java.io.File

import cromwell.binding._
import cromwell.engine.{WorkflowManagerActor, SingleWorkflowRunner}
import cromwell.binding.formatter.{AnsiSyntaxHighlighter, SyntaxFormatter}
import cromwell.parser.WdlParser.SyntaxError
import cromwell.server.CromwellServer
import spray.json._

import scala.util.{Failure, Success}

object Actions extends Enumeration {
  val Parse, Validate, Highlight, Run, Inputs, Server = Value
}

object Main extends App {

  getAction(args.headOption map { _.capitalize }) match {
    case Some(x) if x == Actions.Validate => validate(args.tail)
    case Some(x) if x == Actions.Highlight => highlight(args.tail)
    case Some(x) if x == Actions.Inputs => inputs(args.tail)
    case Some(x) if x == Actions.Run => run(args.tail)
    case Some(x) if x == Actions.Parse => parse(args.tail)
    case Some(x) if x == Actions.Server => CromwellServer
    case None => CromwellServer
    case _ => usageAndExit()
  }

  def validate(args: Array[String]): Unit = {
    if (args.length != 1) usageAndExit()
    try {
      WdlBinding.process(new File(args(0)))
    } catch {
      case e:SyntaxError => println(e)
    }
  }

  def highlight(args: Array[String]) {
    val binding = WdlBinding.process(new File(args(0)))
    val formatter = new SyntaxFormatter(AnsiSyntaxHighlighter)
    println(formatter.format(binding))
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
        runner.workflowManagerActor ! WorkflowManagerActor.Shutdown
      case Failure(f) => println(f.printStackTrace())
    }
  }

  def parse(args: Array[String]): Unit = {
    if (args.length != 1) usageAndExit()
    else println(WdlBinding.getAst(new File(args(0))).toPrettyString)
  }

  def usageAndExit(exit: Boolean = true): Unit = {
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
        |
        |If no action is specified on the command line, the default
        |action is 'server'.
      """.stripMargin)
    if(exit) System.exit(-1)
  }

  def getAction(firstArg: Option[String]): Option[Actions.Value] = for {
    arg <- firstArg
    a <- Actions.values find { _.toString == arg }
  } yield a
}
