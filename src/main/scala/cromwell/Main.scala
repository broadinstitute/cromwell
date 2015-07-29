package cromwell

import java.io.File
import java.nio.file.Paths

import com.typesafe.config.ConfigFactory
import cromwell.binding._
import cromwell.util.FileUtil.{EnhancedFile, EnhancedPath}
import cromwell.binding.formatter.{AnsiSyntaxHighlighter, SyntaxFormatter}
import cromwell.binding.{AstTools, _}
import cromwell.engine.workflow.{WorkflowManagerActor, SingleWorkflowRunnerActor}
import cromwell.parser.WdlParser.SyntaxError
import cromwell.server.{CromwellServer, DefaultWorkflowManagerSystem, WorkflowManagerSystem}
import cromwell.util.FileUtil
import org.slf4j.LoggerFactory
import spray.json._

object Actions extends Enumeration {
  val Parse, Validate, Highlight, Run, Inputs, Server = Value
}

object Main extends App {
  val Properties = System.getProperties
  val LoggerProperty = "CROMWELL_LOGGER"
  lazy val Log = LoggerFactory.getLogger("main")

  Option(Properties.getProperty(LoggerProperty)) match {
    case None => args.headOption.map {_.capitalize}.find {_ == "SERVER"} match {
      case Some(x) => Properties.setProperty(LoggerProperty, "SERVER")
      case _ =>
    }
    case _ =>
  }

  getAction(args.headOption map { _.capitalize }) match {
    case Some(x) if x == Actions.Validate => validate(args.tail)
    case Some(x) if x == Actions.Highlight => highlight(args.tail)
    case Some(x) if x == Actions.Inputs => inputs(args.tail)
    case Some(x) if x == Actions.Run => run(args.tail, DefaultWorkflowManagerSystem())
    case Some(x) if x == Actions.Parse => parse(args.tail)
    case Some(x) if x == Actions.Server => CromwellServer
    case None => usageAndExit(true)
    case _ => usageAndExit()
  }

  def validate(args: Array[String]): Unit = {
    if (args.length != 1) usageAndExit()
    try {
      WdlNamespace.load(new File(args(0)), WorkflowManagerActor.BackendType)
    } catch {
      case e: SyntaxError => println(e)
    }
  }

  def highlight(args: Array[String]) {
    val namespace = WdlNamespace.load(new File(args(0)), WorkflowManagerActor.BackendType)
    val formatter = new SyntaxFormatter(AnsiSyntaxHighlighter)
    println(formatter.format(namespace))
  }

  def inputs(args: Array[String]): Unit = {
    if (args.length != 1) usageAndExit()
    try {
      import cromwell.binding.types.WdlTypeJsonFormatter._
      val namespace = WdlNamespace.load(new File(args(0)), WorkflowManagerActor.BackendType)
      namespace match {
        case x: NamespaceWithWorkflow => println(x.workflow.inputs.toJson.prettyPrint)
        case _ => println("WDL does not have a local workflow")
      }
    } catch {
      case e:SyntaxError => println(e)
    }
  }

  def run(args: Array[String], workflowManagerSystem: WorkflowManagerSystem): Unit = {
    if (args.length != 2) usageAndExit()

    Log.info(s"Backend is: ${WorkflowManagerActor.BackendType}")

    Log.info(s"RUN sub-command")
    Log.info(s"  WDL file: ${args(0)}")
    Log.info(s"  Inputs: ${args(1)}")

    try {
      val wdlSource = Paths.get(args(0)).slurp
      val wdlJson = Paths.get(args(1)).slurp
      val jsValue = wdlJson.parseJson

      val inputs: binding.WorkflowRawInputs = jsValue match {
        case JsObject(rawInputs) => rawInputs
        case _ => throw new RuntimeException("Expecting a JSON object")
      }

      inputs foreach { case (k, v) => Log.info(s"input: $k => $v") }
      val singleWorkflowRunner = SingleWorkflowRunnerActor.props(wdlSource, wdlJson, inputs, workflowManagerSystem.workflowManagerActor)
      workflowManagerSystem.actorSystem.actorOf(singleWorkflowRunner)
      workflowManagerSystem.actorSystem.awaitTermination()
      // And now we just wait for the magic to happen
    } catch {
      case e: Exception =>
        println("Unable to process inputs")
        e.printStackTrace()
    } finally {
      workflowManagerSystem.shutdown() // Future[Unit], but give it a shot
    }
  }

  def parse(args: Array[String]): Unit = {
    if (args.length != 1) usageAndExit()
    else println(AstTools.getAst(new File(args(0))).toPrettyString)
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
      """.stripMargin)
    if (exit) System.exit(-1)
  }

  def getAction(firstArg: Option[String]): Option[Actions.Value] = for {
    arg <- firstArg
    a <- Actions.values find { _.toString == arg }
  } yield a
}
