package cromwell

import java.io.File
import java.nio.file.{Files, Paths}

import cromwell.binding.formatter.{AnsiSyntaxHighlighter, HtmlSyntaxHighlighter, SyntaxFormatter}
import cromwell.binding.{AstTools, _}
import cromwell.engine.WorkflowSourceFiles
import cromwell.engine.workflow.{SingleWorkflowRunnerActor, WorkflowManagerActor, WorkflowOptions}
import cromwell.parser.WdlParser.SyntaxError
import cromwell.server.{CromwellServer, DefaultWorkflowManagerSystem, WorkflowManagerSystem}
import cromwell.util.FileUtil.EnhancedPath
import org.slf4j.LoggerFactory
import spray.json._

import scala.util.{Failure, Success, Try}

object Actions extends Enumeration {
  val Parse, Validate, Highlight, Run, Inputs, Server = Value
}

object Main extends App {
  val Properties = System.getProperties
  val LoggerProperty = "CROMWELL_LOGGER"
  lazy val Log = LoggerFactory.getLogger("main")

  Option(Properties.getProperty(LoggerProperty)) match {
    case None if args.headOption.map {_.toUpperCase}.contains("SERVER") =>
      Properties.setProperty(LoggerProperty, "SERVER")
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
    if (args.length >= 2 && !Seq("html", "console").contains(args(1))) {
      println("usage: highlight <wdl> <html|console>")
      System.exit(1)
    }

    val visitor = args match {
      case a if a.length >= 2 && args(1) == "html" => HtmlSyntaxHighlighter
      case _ => AnsiSyntaxHighlighter
    }
    val namespace = WdlNamespace.load(new File(args(0)), WorkflowManagerActor.BackendType)
    val formatter = new SyntaxFormatter(visitor)
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
    if (args.length < 1 || args.length > 3) usageAndExit()

    args foreach { arg =>
      val path = Paths.get(arg)
      if (arg != "-") {
        if (!Files.exists(path)) {
          System.err.println(s"ERROR: file does not exist: $arg")
          System.exit(1)
        }
        if (!Files.isReadable(path)) {
          System.err.println(s"ERROR: file is not readable: $arg")
          System.exit(1)
        }
      }
    }

    val wdlFile = args(0)
    val inputsJsonFile = Try(args(1)).getOrElse(wdlFile.replaceAll("\\.wdl$", ".json"))
    val workflowOptionsFile = Try(args(2)).getOrElse(wdlFile.replace("\\.wdl$", ".options.json"))

    Log.info(s"Default backend: ${WorkflowManagerActor.BackendType}")
    Log.info(s"RUN sub-command")
    Log.info(s"  WDL file: ${args(0)}")

    try {
      val wdlSource = Paths.get(args(0)).slurp
      val inputsJson = inputsJsonFile match {
        case "-" => "{}"
        case path if path != args(0) && Files.exists(Paths.get(path)) =>
          Log.info(s"  Inputs: $path")
          Paths.get(path).slurp
        case _ =>
          System.err.println(s"ERROR: No workflow inputs specified")
          System.exit(1)
          ""
      }
      val workflowOptions = workflowOptionsFile match {
        case "-" => "{}"
        case path if path != args(0) && Files.exists(Paths.get(path)) =>
          Log.info(s"  Workflow Options: $path")
          Paths.get(path).slurp
        case _ => "{}"
      }
      val jsValue = inputsJson.parseJson

      val inputs: binding.WorkflowRawInputs = jsValue match {
        case JsObject(rawInputs) => rawInputs
        case _ => throw new RuntimeException("Expecting a JSON object")
      }

      val options = WorkflowOptions.fromJsonObject(workflowOptions.parseJson.asInstanceOf[JsObject]) match {
        case Success(x) => x
        case Failure(ex) =>
          System.err.println(s"ERROR: Could not parse workflow options:\n")
          System.err.println(ex.getMessage)
          System.exit(1)
          ???
      }

      inputs foreach { case (k, v) => Log.info(s"input: $k => $v") }
      val sources = WorkflowSourceFiles(wdlSource, inputsJson, options.asPrettyJson)
      val singleWorkflowRunner = SingleWorkflowRunnerActor.props(sources, inputs, workflowManagerSystem.workflowManagerActor)
      workflowManagerSystem.actorSystem.actorOf(singleWorkflowRunner, "SingleWorkflowRunnerActor")
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
        |run <WDL file> [<JSON inputs file> [<JSON workflow options]]
        |
        |  Given a WDL file and JSON file containing the value of the
        |  workflow inputs, this will run the workflow locally and
        |  print out the outputs in JSON format.  The workflow
        |  options file specifies some runtime configuration for the
        |  workflow (see README for details)
        |
        |parse <WDL file>
        |
        |  Compares a WDL file against the grammar and prints out an
        |  abstract syntax tree if it is valid, and a syntax error
        |  otherwise.  Note that higher-level AST checks are not done
        |  via this sub-command and the 'validate' subcommand should
        |  be used for full validation
        |
        |highlight <WDL file> <html|console>
        |
        |  Reformats and colorizes/tags a WDL file. The second
        |  parameter is the output type.  "html" will output the WDL
        |  file with <span> tags around elements.  "console" mode
        |  will output colorized text to the terminal
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
