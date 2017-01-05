package cromwell

import com.typesafe.config.ConfigFactory
import cromwell.engine.workflow.SingleWorkflowRunnerActor
import cromwell.engine.workflow.SingleWorkflowRunnerActor.RunWorkflow
import cromwell.server.{CromwellServer, CromwellSystem}
import cromwell.util.PromiseActor
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.concurrent.{Await, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object Main extends App {
  val CommandLine = CromwellCommandLine(args)
  initLogging(CommandLine)

  lazy val Log = LoggerFactory.getLogger("cromwell")
  lazy val CromwellSystem: CromwellSystem = Try {
    new CromwellSystem {}
  } recoverWith {
    case t: Throwable =>
      Log.error("Failed to instantiate Cromwell System. Shutting down Cromwell.")
      Log.error(t.getMessage)
      System.exit(1)
      Failure(t)
  } get

  CommandLine match {
    case UsageAndExit => usageAndExit()
    case VersionAndExit => versionAndExit()
    case RunServer => waitAndExit(CromwellServer.run(CromwellSystem), CromwellSystem)
    case r: RunSingle => runWorkflow(r)
  }

  /**
    * If a cromwell server is going to be run, makes adjustments to the default logback configuration.
    * Overwrites LOG_MODE system property used in our logback.xml, _before_ the logback classes load.
    * Restored from similar functionality in
    *   https://github.com/broadinstitute/cromwell/commit/2e3f45b#diff-facc2160a82442932c41026c9a1e4b2bL28
    * TODO: Logback is configurable programmatically. We don't have to overwrite system properties like this.
    *
    * Also copies variables from config/system/environment/defaults over to the system properties.
    * Fixes issue where users are trying to specify Java properties as environment variables.
    */
  private def initLogging(commandLine: CromwellCommandLine): Unit = {
    val defaultLogMode = commandLine match {
      case RunServer => "STANDARD"
      case _ => "PRETTY"
    }

    val defaultProps = Map("LOG_MODE" -> defaultLogMode, "LOG_LEVEL" -> "INFO")

    val config = ConfigFactory.load
      .withFallback(ConfigFactory.systemEnvironment())
      .withFallback(ConfigFactory.parseMap(defaultProps.asJava, "Defaults"))

    val props = sys.props
    defaultProps.keys foreach { key =>
      props += key -> config.getString(key)
    }

    /*
    We've possibly copied values from the environment, or our defaults, into the system properties.
    Make sure that the next time one uses the ConfigFactory that our updated system properties are loaded.
     */
    ConfigFactory.invalidateCaches()
  }

  private def runWorkflow(commandLine: RunSingle): Unit = {
    implicit val actorSystem = CromwellSystem.actorSystem

    Log.info(s"RUN sub-command")
    Log.info(s"  WDL file: ${commandLine.wdlPath}")
    commandLine.inputsPath foreach { i => Log.info(s"  Inputs: $i") }
    commandLine.optionsPath foreach { o => Log.info(s"  Workflow Options: $o") }
    commandLine.metadataPath foreach { m => Log.info(s"  Workflow Metadata Output: $m") }

    val runnerProps = SingleWorkflowRunnerActor.props(commandLine.sourceFiles, commandLine.metadataPath)

    val runner = CromwellSystem.actorSystem.actorOf(runnerProps, "SingleWorkflowRunnerActor")

    import PromiseActor.EnhancedActorRef

    waitAndExit(runner.askNoTimeout(RunWorkflow), CromwellSystem)
  }

  private def waitAndExit(futureResult: Future[Any], workflowManagerSystem: CromwellSystem): Unit = {
    Await.ready(futureResult, Duration.Inf)

    workflowManagerSystem.shutdownActorSystem()

    val returnCode = futureResult.value.get match {
      case Success(_) => 0
      case Failure(e) =>
        Console.err.println(e.getMessage)
        1
    }

    sys.exit(returnCode)
  }

  def usageAndExit(): Unit = {
    println(
      """
        |java -jar cromwell.jar <action> <parameters>
        |
        |Actions:
        |run <WDL file> [<JSON inputs file>] [<JSON workflow options>]
        |  [<OUTPUT workflow metadata>] [<Zip of WDL Files>]
        |
        |  Given a WDL file and JSON file containing the value of the
        |  workflow inputs, this will run the workflow locally and
        |  print out the outputs in JSON format.  The workflow
        |  options file specifies some runtime configuration for the
        |  workflow (see README for details).  The workflow metadata
        |  output is an optional file path to output the metadata. The
        |  directory of WDL files is optional. However, it is required
        |  if the primary workflow imports workflows that are outside
        |  of the root directory of the Cromwell project.
        |
        |  Use a single dash ("-") to skip optional files. Ex:
        |    run noinputs.wdl - - metadata.json -
        |
        |server
        |
        |  Starts a web server on port 8000.  See the web server
        |  documentation for more details about the API endpoints.
        |
        |-version
        |
        |   Returns the version of the Cromwell engine.
        |
      """.stripMargin)

    System.exit(1)
  }

  def versionAndExit(): Unit = {
    val versionConf = ConfigFactory.load("cromwell-version.conf").getConfig("version")
    println(
      s"""
         |cromwell: ${versionConf.getString("cromwell")}
       """.stripMargin
    )
    System.exit(1)
  }
}
