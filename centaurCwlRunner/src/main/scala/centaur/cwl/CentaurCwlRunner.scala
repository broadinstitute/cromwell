package centaur.cwl

import better.files._
import centaur.api.CentaurCromwellClient
import centaur.cwl.Outputs._
import centaur.test.TestOptions
import centaur.test.standard.{CentaurTestCase, CentaurTestFormat}
import centaur.test.submit.{SubmitHttpResponse, SubmitWorkflowResponse}
import centaur.test.workflow.{AllBackendsRequired, Workflow, WorkflowData}
import com.google.api.services.storage.StorageScopes
import com.typesafe.scalalogging.StrictLogging
import common.util.VersionUtil
import cromwell.api.model.{Aborted, Failed, NonTerminalStatus, Succeeded}
import cromwell.cloudsupport.gcp.auth.{ApplicationDefaultMode, ServiceAccountMode}
import cromwell.cloudsupport.gcp.auth.ServiceAccountMode.JsonFileFormat
import cromwell.core.WorkflowOptions
import cromwell.core.path.{DefaultPathBuilderFactory, PathBuilderFactory}
import cromwell.filesystems.gcs.GcsPathBuilderFactory
import spray.json._

import scala.concurrent.Await
import scala.concurrent.duration.Duration

/**
  * Runs workflows in a "cwl-runner" friendly way.
  *
  * https://github.com/broadinstitute/cromwell/issues/2590
  * https://github.com/common-workflow-language/common-workflow-language/blob/v1.0.1/CONFORMANCE_TESTS.md
  * https://github.com/common-workflow-language/common-workflow-language/blob/v1.0.1/draft-3/cwl-runner.cwl#L5-L68
  * https://github.com/common-workflow-language/common-workflow-language/pull/278/files#diff-ee814a9c027fc9750beb075c283a973cR49
  */
object CentaurCwlRunner extends StrictLogging {

  case class CommandLineArguments(workflowSource: Option[File] = None,
                                  workflowInputs: Option[File] = None,
                                  quiet: Boolean = false,
                                  outdir: Option[File] = None,
                                  localMode: Boolean = true,
                                  serviceAccountJson: Option[String] = None)

  // TODO: This would be cleaner with Enumeratum
  object ExitCode extends Enumeration {

    protected case class Val(status: Int) extends super.Val

    implicit class ValueToVal(val exitCodeValue: Value) extends AnyVal {
      def status: Int = exitCodeValue.asInstanceOf[Val].status
    }

    val Success = Val(0)
    val Failure = Val(1)
    val NotImplemented = Val(33)
  }

  private val parser = buildParser()
  private lazy val centaurCwlRunnerVersion = VersionUtil.getVersion("centaur-cwl-runner")

  private def showUsage(): ExitCode.Value = {
    parser.showUsage()
    ExitCode.Failure
  }

  private def buildParser(): scopt.OptionParser[CommandLineArguments] = {
    new scopt.OptionParser[CommandLineArguments]("java -jar /path/to/centaurCwlRunner.jar") {
      head("centaur-cwl-runner", centaurCwlRunnerVersion)

      help("help").text("Centaur CWL Runner - Cromwell integration testing environment")

      version("version").text("Print version and exit")

      arg[String]("workflow-source").text("Workflow source file.").required().
        action((s, c) => c.copy(workflowSource = Option(File(s))))

      arg[String]("inputs").text("Workflow inputs file.").optional().
        action((s, c) => c.copy(workflowInputs = Option(File(s))))

      opt[Unit]("quiet").text("Only print warnings and errors.").optional().
        action((_, c) => c.copy(quiet = true))

      opt[String]("outdir").text("Output directory, default current directory.").optional().
        action((s, c) =>
          c.copy(outdir = Option(File(s))))

      opt[Unit]("mode-local").text("Assume local file paths (default)").optional().
        action((_, c) => c.copy(localMode = true))

      opt[Unit]("mode-gcs-default").text("Assume GCS file paths, use application default credentials").optional().
        action((_, c) => c.copy(localMode = false))

      opt[String]("mode-gcs-service").text("Assume GCS file paths. Use service account credentials, path to credentials JSON as an argument.").
        optional().
        action((s, c) => c.copy(localMode = false, serviceAccountJson = Option(s)))
    }
  }

  private def runCentaur(args: CommandLineArguments): ExitCode.Value = {

    def zipSiblings(file: File): File = {
      val zipFile = File.newTemporaryFile("cwl_imports.", ".zip")
      val dir = file.parent
      if (!args.quiet) {
        logger.info(s"Zipping $dir to $zipFile")
      }
      val files = dir.children
      Cmds.zip(files.toSeq: _*)(zipFile)
      zipFile
    }

    val workflowPath = args.workflowSource.get
    val (parsedWorkflowPath, workflowRoot) = workflowPath.path.toAbsolutePath.toString.split("#") match {
      case Array(file) => File(file) -> None
      case Array(file, root) => File(file) -> Option(root)
    }
    val outdirOption = args.outdir.map(_.pathAsString)
    val testName = workflowPath.name
    val workflowContents = parsedWorkflowPath.contentAsString
    val inputContents = args.workflowInputs.map(_.contentAsString)
    val workflowType = Option("cwl")
    val workflowTypeVersion = None
    val optionsContents = outdirOption map { outdir =>
      JsObject("cwl_outdir" -> JsString(outdir)).compactPrint
    }
    val labels = List.empty
    val zippedImports = Option(zipSiblings(workflowPath)) // TODO: Zipping all the things! Be more selective.
    val backends = AllBackendsRequired(List.empty)
    val workflowMetadata = None
    val notInMetadata = List.empty
    val directoryContentCounts = None
    val testFormat = CentaurTestFormat.WorkflowSuccessTest
    val testOptions = TestOptions(List.empty, ignore = false)
    val submitResponseOption = None

    val workflowData = WorkflowData(
      workflowContents, workflowRoot, workflowType, workflowTypeVersion, inputContents, optionsContents, labels, zippedImports)
    val workflow = Workflow(testName, workflowData, workflowMetadata, notInMetadata, directoryContentCounts, backends)
    val testCase = CentaurTestCase(workflow, testFormat, testOptions, submitResponseOption)

    if (!args.quiet) {
      logger.info(s"Starting test for $workflowPath")
    }

    def pathBuilderFactory: PathBuilderFactory = {
      (args.localMode, args.serviceAccountJson) match {
        case (true, _) => DefaultPathBuilderFactory
        case (false, Some(serviceAccountJson)) =>
          import scala.collection.JavaConverters._
          val scopes = List(
            StorageScopes.DEVSTORAGE_FULL_CONTROL,
            StorageScopes.DEVSTORAGE_READ_WRITE
          ).asJava
          val auth = ServiceAccountMode("service-account", JsonFileFormat(serviceAccountJson), scopes)
          GcsPathBuilderFactory(auth, "centaur-cwl", None)
        case (false, _) =>
          GcsPathBuilderFactory(ApplicationDefaultMode("application-default"), "centaur-cwl", None)
      }
    }

    try {
      import CentaurCromwellClient.{blockingEc, system}
      lazy val pathBuilder = Await.result(pathBuilderFactory.withOptions(WorkflowOptions.empty), Duration.Inf)

      testCase.testFunction.run.get match {
        case unexpected: SubmitHttpResponse =>
          logger.error(s"Unexpected response: $unexpected")
          ExitCode.Failure
        case SubmitWorkflowResponse(submittedWorkflow) =>
          val status = CentaurCromwellClient.status(submittedWorkflow).get
          status match {
            case unexpected: NonTerminalStatus =>
              logger.error(s"Unexpected status: $unexpected")
              ExitCode.Failure
            case Aborted =>
              logger.error(s"Unexpected abort.")
              ExitCode.Failure
            case Failed =>
              logger.error(s"Unexpected failure.")
              ExitCode.Failure
            case Succeeded =>
              val result = handleOutput(submittedWorkflow, pathBuilder)
              if (!args.quiet) {
                logger.info(s"Result: $result")
              } else {
                // Print directly to stdout during --quiet
                println(result)
              }
              ExitCode.Success
          }
      }
    } finally {
      zippedImports map { zipFile =>
        if (!args.quiet) {
          logger.info(s"Deleting $zipFile")
        }
        zipFile.delete(swallowIOExceptions = true)
      }
      Await.result(CentaurCromwellClient.system.terminate(), Duration.Inf)
      ()
    }
  }

  def main(args: Array[String]): Unit = {
    val parsedArgsOption = parser.parse(args, CommandLineArguments())
    val exitCode = parsedArgsOption match {
      case Some(parsedArgs) => runCentaur(parsedArgs)
      case None => showUsage()
    }
    System.exit(exitCode.status)
  }
}
