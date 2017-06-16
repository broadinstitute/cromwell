package cromwell

import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._
import com.typesafe.config.ConfigFactory
import cromwell.CommandLineParser._
import cromwell.core.path.Path
import cromwell.core.{WorkflowSourceFilesCollection, WorkflowSourceFilesWithDependenciesZip, WorkflowSourceFilesWithoutImports}
import cromwell.engine.workflow.SingleWorkflowRunnerActor
import cromwell.engine.workflow.SingleWorkflowRunnerActor.RunWorkflow
import cromwell.server.{CromwellServer, CromwellSystem}
import lenthall.exception.MessageAggregation
import lenthall.validation.ErrorOr._
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.concurrent.duration.{Duration, _}
import scala.concurrent.{Await, Future, TimeoutException}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}


object CromwellEntryPoint {

  /**
    * Run Cromwell in server mode.
    */
  def runServer() = {
    val system = buildCromwellSystem(Server)
    waitAndExit(CromwellServer.run, system)
  }

  /**
    * Run a single workflow using the successfully parsed but as yet not validated arguments.
    */
  def runSingle(args: CommandLineArguments): Unit = {
    val cromwellSystem = buildCromwellSystem(Run)
    implicit val actorSystem = cromwellSystem.actorSystem

    val sources = validateRunArguments(args)
    val runnerProps = SingleWorkflowRunnerActor.props(sources, args.metadataOutput)(cromwellSystem.materializer)

    val runner = cromwellSystem.actorSystem.actorOf(runnerProps, "SingleWorkflowRunnerActor")

    import cromwell.util.PromiseActor.EnhancedActorRef
    waitAndExit(_ => runner.askNoTimeout(RunWorkflow), cromwellSystem)
  }

  private def buildCromwellSystem(command: Command): CromwellSystem = {
    initLogging(command)
    lazy val Log = LoggerFactory.getLogger("cromwell")
    Try {
      new CromwellSystem {}
    } recoverWith {
      case t: Throwable =>
        Log.error("Failed to instantiate Cromwell System. Shutting down Cromwell.")
        Log.error(t.getMessage)
        System.exit(1)
        Failure(t)
    } get
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
  private def initLogging(command: Command): Unit = {
    val logbackSetting = command match {
      case Server => "STANDARD"
      case Run => "PRETTY"
    }

    val defaultProps = Map(
      "LOG_MODE" -> logbackSetting,
      "LOG_LEVEL" -> "INFO"
    )

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

  private def waitAndExit(runner: CromwellSystem => Future[Any], workflowManagerSystem: CromwellSystem): Unit = {
    val futureResult = runner(workflowManagerSystem)
    Await.ready(futureResult, Duration.Inf)

    try {
      Await.ready(workflowManagerSystem.shutdownActorSystem(), 30 seconds)
    } catch {
      case _: TimeoutException => Console.err.println("Timed out trying to shutdown actor system")
      case other: Exception => Console.err.println(s"Unexpected error trying to shutdown actor system: ${other.getMessage}")
    }

    val returnCode = futureResult.value.get match {
      case Success(_) => 0
      case Failure(e) =>
        Console.err.println(e.getMessage)
        1
    }

    sys.exit(returnCode)
  }

  def validateRunArguments(args: CommandLineArguments): WorkflowSourceFilesCollection = {

    val workflowSource = readContent("Workflow source", args.workflowSource.get)
    val inputsJson = readJson("Workflow inputs", args.workflowInputs)
    val optionsJson = readJson("Workflow options", args.workflowOptions)
    val labelsJson = readJson("Workflow labels", args.workflowLabels)

    val sourceFileCollection = args.imports match {
      case Some(p) => (workflowSource |@| inputsJson |@| optionsJson |@| labelsJson) map { (w, i, o, l) =>
        WorkflowSourceFilesWithDependenciesZip.apply(
          workflowSource = w,
          workflowType = Option("WDL"),
          workflowTypeVersion = None,
          inputsJson = i,
          workflowOptionsJson = o,
          labelsJson = l,
          importsZip = p.loadBytes)
      }
      case None => (workflowSource |@| inputsJson |@| optionsJson |@| labelsJson) map { (w, i, o, l) =>
        WorkflowSourceFilesWithoutImports.apply(
          workflowSource = w,
          workflowType = Option("WDL"),
          workflowTypeVersion = None,
          inputsJson = i,
          workflowOptionsJson = o,
          labelsJson = l
        )
      }
    }

    val sourceFiles = for {
      sources <- sourceFileCollection
      _ <- writeableMetadataPath(args.metadataOutput)
    } yield sources

    sourceFiles match {
      case Valid(r) => r
      case Invalid(nel) => throw new RuntimeException with MessageAggregation {
        override def exceptionContext: String = "ERROR: Unable to run Cromwell:"
        override def errorMessages: Traversable[String] = nel.toList
      }
    }
  }

  private def writeableMetadataPath(path: Option[Path]): ErrorOr[Unit] = {
    path match {
      case Some(p) if !metadataPathIsWriteable(p) => s"Unable to write to metadata directory: $p".invalidNel
      case _ => ().validNel
    }
  }

  /** Read the path to a string. */
  private def readContent(inputDescription: String, path: Path): ErrorOr[String] = {
    if (!path.exists) {
      s"$inputDescription does not exist: $path".invalidNel
    } else if (!path.isReadable) {
      s"$inputDescription is not readable: $path".invalidNel
    } else path.contentAsString.validNel
  }

  /** Read the path to a string, unless the path is None, in which case returns "{}". */
  private def readJson(inputDescription: String, pathOption: Option[Path]): ErrorOr[String] = {
    pathOption match {
      case Some(path) => readContent(inputDescription, path)
      case None => "{}".validNel
    }
  }

  private def metadataPathIsWriteable(metadataPath: Path): Boolean =
    Try(metadataPath.createIfNotExists(createParents = true).append("")).isSuccess

}
