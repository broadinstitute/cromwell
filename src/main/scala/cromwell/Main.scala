package cromwell

import java.io.{File => JFile}
import java.nio.file.{Files, Path, Paths}

import akka.actor.{Actor, Props, Status}
import better.files._
import cromwell.util.TerminalUtil
import wdl4s.formatter.{AnsiSyntaxHighlighter, HtmlSyntaxHighlighter, SyntaxFormatter}
import wdl4s.{AstTools, _}
import cromwell.engine.WorkflowSourceFiles
import cromwell.engine.workflow.SingleWorkflowRunnerActor.RunWorkflow
import cromwell.engine.workflow.{SingleWorkflowRunnerActor, WorkflowOptions}
import cromwell.instrumentation.Instrumentation.Monitor
import cromwell.server.{CromwellServer, WorkflowManagerSystem}
import cromwell.util.FileUtil._
import org.slf4j.LoggerFactory
import spray.json._

import scala.concurrent.duration._
import scala.concurrent.{Future, Await, Promise}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import Main.deprecationNotice

object Actions extends Enumeration {
  val Parse, Validate, Highlight, Run, Inputs, Server = Value
}

object Main extends App {
  setupServerLogging(args)

  /*
   * scala.App's DelayedInit is tricky, as the docs say. During tests we definitely don't want to use sys.exit on an
   * error, and while testing "run ..." we want to change to a test workflow manager system. Unfortunately...
   *
   * - http://stackoverflow.com/questions/17716975/weird-initialisation-when-override-val-met-delayedinit-in-scala
   * - http://stackoverflow.com/questions/13009265/what-happens-in-scala-when-loading-objects-that-extends-app
   * - http://stackoverflow.com/questions/18582709/so-what-can-you-use-delayedinit-for-safely
   * - http://stackoverflow.com/questions/24437423/in-scala-should-i-use-the-app-trait
   * - etc. etc. etc.
   *
   * So, split the main logic out so that it's easier to change these values in tests.
   *
   * Also now passing args to runAction instead of the constructor, as even sbt seemed to have issues with the args
   * array becoming null in "new Main(args)" when used with: sbt 'run run ...'
   */
  sys.exit(new Main().runAction(args))

  /**
    * If a cromwell server is going to be run, makes adjustments to the default logback configuration.
    * Overwrites LOG_MODE system property used in our logback.xml, _before_ the logback classes load.
    * Restored from similar functionality in
    *   https://github.com/broadinstitute/cromwell/commit/2e3f45b#diff-facc2160a82442932c41026c9a1e4b2bL28
    * TODO: Logback is configurable programmatically. We don't have to overwrite system properties like this.
    *
    * @param args The command line arguments.
    */
  private def setupServerLogging(args: Array[String]): Unit = {
    getAction(args) match {
      case Some(Actions.Server) => sys.props.getOrElseUpdate("LOG_MODE", "STANDARD")
      case _ =>
    }
  }

  private def getAction(args: Seq[String]): Option[Actions.Value] = for {
    arg <- args.headOption
    argCapitalized = arg.capitalize
    action <- Actions.values find (_.toString == argCapitalized)
  } yield action

  def deprecationNotice(): Unit = {
    Console.err.println(TerminalUtil.highlight(1, "This functionality is deprecated and will be removed in 0.18. Please use wdltool: https://github.com/broadinstitute/wdltool"))
  }
}

/** A simplified version of the Akka `PromiseActorRef` that doesn't time out. */
private class PromiseActor(promise: Promise[Any]) extends Actor {
  override def receive = {
    case Status.Failure(f) =>
      promise.tryFailure(f)
      context.stop(self)
    case success =>
      promise.trySuccess(success)
      context.stop(self)
  }
}

class Main private[cromwell](managerSystem: WorkflowManagerSystem) {
  private[cromwell] val initLoggingReturnCode = initLogging()

  lazy val Log = LoggerFactory.getLogger("cromwell")
  Monitor.start()

  def this() = this(managerSystem = new WorkflowManagerSystem {})

  // CromwellServer still doesn't clean up... so => Any
  def runAction(args: Seq[String]): Int = {
    Main.getAction(args) match {
      case Some(x) if x == Actions.Validate => validate(args.tail)
      case Some(x) if x == Actions.Highlight => highlight(args.tail)
      case Some(x) if x == Actions.Inputs => inputs(args.tail)
      case Some(x) if x == Actions.Run => run(args.tail)
      case Some(x) if x == Actions.Parse => parse(args.tail)
      case Some(x) if x == Actions.Server => runServer(args.tail)
      case _ => usageAndExit()
    }
  }

  def validate(args: Seq[String]): Int = {
    deprecationNotice()
    continueIf(args.length == 1) {
      loadWdl(args.head) { _ => 0 }
    }
  }

  def highlight(args: Seq[String]): Int = {
    deprecationNotice()
    continueIf(args.length == 2 && Seq("html", "console").contains(args(1))) {
      loadWdl(args.head) { namespace =>
        val formatter = new SyntaxFormatter(if (args(1) == "html") HtmlSyntaxHighlighter else AnsiSyntaxHighlighter)
        println(formatter.format(namespace))
        0
      }
    }
  }

  def inputs(args: Seq[String]): Int = {
    deprecationNotice()
    continueIf(args.length == 1) {
      loadWdl(args.head) { namespace =>
        import wdl4s.types.WdlTypeJsonFormatter._
        namespace match {
          case x: NamespaceWithWorkflow => println(x.workflow.inputs.toJson.prettyPrint)
          case _ => println("WDL does not have a local workflow")
        }
        0
      }
    }
  }

  def runServer(args: Seq[String]): Int = {
    continueIf(args.isEmpty)(waitAndExit(CromwellServer.run(), CromwellServer))
  }

  /* Begin .run() method and utilities */

  val WdlLabel = "WDL file"
  val InputsLabel = "Inputs"
  val OptionsLabel = "Workflow Options"
  val MetadataLabel = "Workflow Metadata Output"

  def run(args: Seq[String]): Int = {
    continueIf(args.nonEmpty && args.length <= 4) {
      val wdlPath = Paths.get(args.head)
      val inputsPath = argPath(args, 1, Option(".inputs"), checkDefaultExists = false)
      val optionsPath = argPath(args, 2, Option(".options"), checkDefaultExists = true)
      val metadataPath = argPath(args, 3, None)

      Log.info(s"RUN sub-command")
      Log.info(s"  $WdlLabel: $wdlPath")
      inputsPath.foreach(path => Log.info(s"  $InputsLabel: $path"))
      optionsPath.foreach(path => Log.info(s"  $OptionsLabel: $path"))
      metadataPath.foreach(path => Log.info(s"  $MetadataLabel: $path"))

      checkPaths(wdlPath, inputsPath, optionsPath, metadataPath) match {
        case Success(workflowSourceFiles) => runWorkflow(workflowSourceFiles, metadataPath)
        case Failure(ex) =>
          Console.err.println(ex.getMessage)
          1
      }
    }
  }

  /**
    * Retrieve the arg at index as path, or return some default. Args specified as "-" will be returned as None.
    *
    * @param args The run command arguments, with the wdl path at arg.head.
    * @param index The index of the path we're looking for.
    * @param defaultExt The default extension to use if the argument was not specified at all.
    * @param checkDefaultExists If true, verify that our computed default file exists before using it.
    * @return The argument as a Path resolved as a sibling to the wdl path.
    */
  private[this] def argPath(args: Seq[String], index: Int, defaultExt: Option[String],
                            checkDefaultExists: Boolean = true): Option[Path] = {

    // To return a default, swap the extension, and then maybe check if the file exists.
    def defaultPath = defaultExt
      .map(ext => swapExt(args.head, ".wdl", ext))
      .filter(path => !checkDefaultExists || Files.exists(Paths.get(path)))

    // Return the path for the arg index, or the default, but remove "-" paths.
    for {
      path <- args.lift(index) orElse defaultPath filterNot (_ == "-")
    } yield Paths.get(path)
  }

  private[this] def runWorkflow(workflowSourceFiles: WorkflowSourceFiles, metadataPath: Option[Path]): Int = {
    val workflowManagerSystem = managerSystem
    val runnerProps = SingleWorkflowRunnerActor.props(workflowSourceFiles, metadataPath,
      workflowManagerSystem.workflowManagerActor)
    val runner = workflowManagerSystem.actorSystem.actorOf(runnerProps, "SingleWorkflowRunnerActor")

    val promise = Promise[Any]()
    val promiseActor = workflowManagerSystem.actorSystem.actorOf(Props(classOf[PromiseActor], promise))
    runner.tell(RunWorkflow, promiseActor)
    waitAndExit(promise.future, workflowManagerSystem)
  }

  private[this] def waitAndExit(futureResult: Future[Any], workflowManagerSystem: WorkflowManagerSystem): Int = {
    Await.ready(futureResult, Duration.Inf)

    workflowManagerSystem.shutdownActorSystem()

    futureResult.value.get match {
      case Success(_) => 0
      case Failure(e) =>
        Console.err.println(e.getMessage)
        1
    }
  }

  /* Utilities for the mini-DSL used in tryWorkflowSourceFiles. */

  /** Read the path to a string. */
  private[this] def readContent(inputDescription: String, path: Path): String = {
    // Provide slightly better errors than the java.io.IOException messages.
    if (!Files.exists(path)) {
      throw new RuntimeException(s"$inputDescription does not exist: $path")
    }
    if (!Files.isReadable(path)) {
      throw new RuntimeException(s"$inputDescription is not readable: $path")
    }
    path.contentAsString
  }

  /** Read the path to a string, unless the path is None, in which case returns "{}". */
  private[this] def readJson(inputDescription: String, pathOption: Option[Path]): String = {
    pathOption match {
      case Some(path) => readContent(inputDescription, path)
      case None => "{}"
    }
  }

  /** Try to write to the path by appending a blank string to file. */
  private[this] def writeTo(outputDescription: String, pathOption: Option[Path]): Unit = {
    pathOption.foreach(_.createIfNotExists().append(""))
  }

  /**
    * Run sanity checks on each of the argument paths.
    *
    * Check if we can read the wdl, inputs, and options Path => String.
    * Then test if we can parse the json files, especially passing the options json to
    * WorkflowOptions.fromJsonObject for special handling.
    * Lastly, ensure we can write to the metadata file.
    *
    * If all is ok, return a Success(WorkflowSourceFiles), otherwise a Failure describing the error.
    */
  private[this] def checkPaths(wdlPath: Path, inputsPath: Option[Path], optionsPath: Option[Path],
                               metadataPath: Option[Path]): Try[WorkflowSourceFiles] = {
    /*
     * NOTE: Check with others regarding backwards compatibility, then feel free to rip out and replace this DSL if it
     * would be faster or cleaner just to re-implement!
     */
    import MainRunDsl._
    for {
      wdlSource <- Trying.to(readContent _, wdlPath, ShouldProcess, WdlLabel)
      inputsJson <- Trying.to(readJson _, inputsPath, ShouldProcess, InputsLabel)
      workflowOptions <- Trying.to(readJson _, optionsPath, ShouldProcess, OptionsLabel)
      _ <- Trying.to(writeTo _, metadataPath, ShouldAccess, MetadataLabel)
    } yield WorkflowSourceFiles(wdlSource, inputsJson, workflowOptions)
  }

  /**
    * Run the workflow options json string through the WorkflowOptions parser.
    * NOTE: It's possible WorkflowOptions.fromJsonObject is being called twice, first here, and secondly within the
    * WorkflowActor system.
    */
  def parseOptions(workflowOptions: String): Try[String] = {
    WorkflowOptions.fromJsonObject(workflowOptions.parseJson.asInstanceOf[JsObject]).map(_.asPrettyJson)
  }

  def parse(args: Seq[String]): Int = {
    deprecationNotice()
    continueIf(args.length == 1) {
      println(AstTools.getAst(new JFile(args.head)).toPrettyString)
      0
    }
  }

  def usageAndExit(): Int = {
    println(
      """
        |java -jar cromwell.jar <action> <parameters>
        |
        |Actions:
        |run <WDL file> [<JSON inputs file> [<JSON workflow options>
        |  [<OUTPUT workflow metadata>]]]
        |
        |  Given a WDL file and JSON file containing the value of the
        |  workflow inputs, this will run the workflow locally and
        |  print out the outputs in JSON format.  The workflow
        |  options file specifies some runtime configuration for the
        |  workflow (see README for details).  The workflow metadata
        |  output is an optional file path to output the metadata.
        |  Use a single dash ("-") to skip optional files. Ex:
        |    run noinputs.wdl - - metadata.json
        |
        |server
        |
        |  Starts a web server on port 8000.  See the web server
        |  documentation for more details about the API endpoints.
        |
        |parse <WDL file>
        |
        |  This functionality is deprecated. Please use wdltool: https://github.com/broadinstitute/wdltool
        |
        |  Compares a WDL file against the grammar and prints out an
        |  abstract syntax tree if it is valid, and a syntax error
        |  otherwise.  Note that higher-level AST checks are not done
        |  via this sub-command and the 'validate' subcommand should
        |  be used for full validation
        |
        |validate <WDL file>
        |
        |  This functionality is deprecated. Please use wdltool: https://github.com/broadinstitute/wdltool
        |
        |  Performs full validation of the WDL file including syntax
        |  and semantic checking
        |
        |inputs <WDL file>
        |
        |  This functionality is deprecated. Please use wdltool: https://github.com/broadinstitute/wdltool
        |
        |  Print a JSON skeleton file of the inputs needed for this
        |  workflow.  Fill in the values in this JSON document and
        |  pass it in to the 'run' subcommand.
        |
        |highlight <WDL file> <html|console>
        |
        |  This functionality is deprecated. Please use wdltool: https://github.com/broadinstitute/wdltool
        |
        |  Reformats and colorizes/tags a WDL file. The second
        |  parameter is the output type.  "html" will output the WDL
        |  file with <span> tags around elements.  "console" mode
        |  will output colorized text to the terminal
        |
      """.stripMargin)
    -1
  }

  private[this] def initLogging(): Int = {
    val systemProperties = sys.props
    systemProperties.getOrElseUpdate("LOG_MODE", "PRETTY")
    systemProperties.getOrElseUpdate("LOG_LEVEL", "INFO")
    /*
    TODO: When implementing DSDEEPB-2271, make sure not to regress DSDEEPB-2339.
    If a parameter is not going to be used due to some combination of arguments, log a warning.
     */
    if (systemProperties.get("LOG_ROOT").isDefined) {
      Console.err.println("WARNING: The LOG_ROOT parameter is currently ignored.")
    }
    val logRoot = systemProperties.getOrElseUpdate("LOG_ROOT", File(".").fullPath)

    try {
      File(logRoot).createDirectories()
      0
    } catch {
      case e: Throwable =>
        Console.err.println(s"Could not create log directory: $logRoot")
        e.printStackTrace()
        1
    }
  }

  private[this] def continueIf(valid: => Boolean)(block: => Int): Int = if (valid) block else usageAndExit()

  private[this] def loadWdl(path: String)(f: WdlNamespace => Int): Int = {
    Try(WdlNamespace.load(new JFile(path))) match {
      case Success(namespace) => f(namespace)
      case Failure(t) =>
        println(t.getMessage)
        1
    }
  }
}
