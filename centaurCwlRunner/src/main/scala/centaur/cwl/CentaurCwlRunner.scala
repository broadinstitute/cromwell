package centaur.cwl

import better.files._
import centaur.api.CentaurCromwellClient
import centaur.test.TestOptions
import centaur.test.standard.{CentaurTestCase, CentaurTestFormat}
import centaur.test.submit.{SubmitHttpResponse, SubmitWorkflowResponse}
import centaur.test.workflow.{AllBackendsRequired, Workflow, WorkflowData}
import com.typesafe.config.ConfigFactory
import cromwell.api.model.{Aborted, Failed, NonTerminalStatus, Succeeded}
import spray.json._

/**
  * Runs workflows in a "cwl-runner" friendly way.
  *
  * https://github.com/broadinstitute/cromwell/issues/2590
  * https://github.com/common-workflow-language/common-workflow-language/blob/v1.0.1/CONFORMANCE_TESTS.md
  * https://github.com/common-workflow-language/common-workflow-language/blob/v1.0.1/draft-3/cwl-runner.cwl#L5-L68
  * https://github.com/common-workflow-language/common-workflow-language/pull/278/files#diff-ee814a9c027fc9750beb075c283a973cR49
  */
object CentaurCwlRunner {

  case class CommandLineArguments(workflowSource: Option[File] = None,
                                  workflowInputs: Option[File] = None,
                                  quiet: Boolean = false,
                                  outdir: Option[File] = None)

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
  private lazy val centaurVersion = ConfigFactory.load("cromwell-version.conf").getString("version.cromwell")

  private def showUsage(): ExitCode.Value = {
    parser.showUsage()
    ExitCode.Failure
  }

  private def buildParser(): scopt.OptionParser[CommandLineArguments] = {
    new scopt.OptionParser[CommandLineArguments]("java -jar /path/to/centaurCwlRunner.jar") {
      head("centaur-cwl-runner", centaurVersion)

      help("help").text("Centaur CWL Runner - Cromwell integration testing environment")

      version("version").text("Print version and exit")

      arg[String]("workflow-source").text("Workflow source file.").required().
        action((s, c) => c.copy(workflowSource = Option(File(s))))

      arg[String]("inputs").text("Workflow inputs file.").optional().
        action((s, c) => c.copy(workflowInputs = Option(File(s))))

      opt[Unit]("quiet").text("Only print warnings and errors.").optional().
        action((_, c) => c.copy(quiet = true))

      opt[String]("outdir").text("Output directory, default current directory. Currently ignored.").optional().
        action((s, c) =>
          c.copy(outdir = Option(File(s))))
    }
  }

  private def runCentaur(args: CommandLineArguments): ExitCode.Value = {

    def zipSiblings(file: File): File = {
      val zipFile = File.newTemporaryFile("cwl_imports.", ".zip")
      val dir = file.parent
      if (!args.quiet) {
        println(s"Zipping $dir to $zipFile")
      }
      val files = dir.children
      Cmds.zip(files.toSeq: _*)(zipFile)
      zipFile
    }

    val workflowPath = args.workflowSource.get
    val outdirOption = args.outdir.map(_.pathAsString)
    val testName = workflowPath.name
    val workflowContents = workflowPath.contentAsString
    val inputContents = args.workflowInputs.map(_.contentAsString)
    val workflowType = workflowPath.extension(includeDot = false)
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
      workflowContents, workflowType, workflowTypeVersion, inputContents, optionsContents, labels, zippedImports)
    val workflow = Workflow(testName, workflowData, workflowMetadata, notInMetadata, directoryContentCounts, backends)
    val testCase = CentaurTestCase(workflow, testFormat, testOptions, submitResponseOption)

    if (!args.quiet) {
      println(s"Starting test for $workflowPath")
    }

    try {
      testCase.testFunction.run.get match {
        case unexpected: SubmitHttpResponse =>
          println(s"Unexpected response: $unexpected")
          ExitCode.Failure
        case SubmitWorkflowResponse(submittedWorkflow) =>
          val status = CentaurCromwellClient.status(submittedWorkflow).get
          status match {
            case unexpected: NonTerminalStatus =>
              println(s"Unexpected status: $unexpected")
              ExitCode.Failure
            case Aborted =>
              println(s"Unexpected abort.")
              ExitCode.Failure
            case Failed =>
              println(s"Unexpected failure.")
              ExitCode.Failure
            case Succeeded =>
              val outputs = CentaurCromwellClient.outputs(submittedWorkflow).get.outputs
              if (!args.quiet) {
                println("Result:")
              }
              println(outputs)
              ExitCode.Success
          }
      }
    } finally {
      zippedImports map { zipFile =>
        if (!args.quiet) {
          println(s"Deleting $zipFile")
        }
        zipFile.delete(swallowIOExceptions = true)
      }
      CentaurCromwellClient.system.terminate()
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
