package cromwell.backend.impl.local

import java.nio.file.{FileSystems, Path, Paths}

import akka.actor.Props
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionAbortedResponse, BackendJobExecutionFailedResponse, BackendJobExecutionResponse, BackendJobExecutionSucceededResponse}
import cromwell.backend._
import org.slf4j.LoggerFactory
import wdl4s._
import wdl4s.util.TryUtil
import wdl4s.values.{WdlFile, WdlValue}

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object LocalJobExecutionActor {
  val ProcessKilledCode = 143
  val logger = LoggerFactory.getLogger("LocalBackend")
  // TODO Support GCS ?
  val fileSystems = List(FileSystems.getDefault)

  private def splitFqn(fullyQualifiedName: FullyQualifiedName): (String, String) = {
    val lastIndex = fullyQualifiedName.lastIndexOf(".")
    (fullyQualifiedName.substring(0, lastIndex), fullyQualifiedName.substring(lastIndex + 1))
  }

  case class Command(argv: Seq[String]) {
    override def toString = argv.map(s => "\"" + s + "\"").mkString(" ")
  }

  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor): Props =
    Props(new LocalJobExecutionActor(jobDescriptor, configurationDescriptor))
}

class LocalJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                             override val configurationDescriptor: BackendConfigurationDescriptor)
  extends BackendJobExecutionActor with SharedFileSystem {

  import LocalJobExecutionActor._
  import better.files._
  import cromwell.core.PathFactory._

  import scala.sys.process._

  override implicit val ec: ExecutionContext = scala.concurrent.ExecutionContext.global

  // Mutable process variable assigned when the execute method is called.
  // Needs to be accessible to be killed when the job is aborted.
  private var process: Option[Process] = None

  val workflowDescriptor = jobDescriptor.descriptor
  val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobDescriptor.key)
  val fileSystemsConfig = configurationDescriptor.backendConfig.getConfig("filesystems")
  override val sharedFsConfig = fileSystemsConfig.getConfig("local")

  val call = jobDescriptor.key.call
  val callEngineFunction =  LocalJobExpressionFunctions(jobPaths)

  val lookup = jobDescriptor.inputs.apply _

  private def evaluate(wdlExpression: WdlExpression) = wdlExpression.evaluate(lookup, callEngineFunction)

  val runtimeAttributes = {
    val evaluateAttrs = call.task.runtimeAttributes.attrs mapValues evaluate
    // Fail the call if runtime attributes can't be evaluated
    val evaluatedAttributes = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    LocalRuntimeAttributes(evaluatedAttributes)
  }

  lazy val runsOnDocker = runtimeAttributes.dockerImage.isDefined
  lazy val processArgs = {
    val dockerRun = runtimeAttributes.dockerImage map buildDockerRunCommand getOrElse ""
    Command(Seq("/bin/bash", "-c", s"cat ${jobPaths.script} | $dockerRun /bin/bash <&0"))
  }

  // Stream Writers
  lazy val stdoutWriter = jobPaths.stdout.untailed
  lazy val stderrTailed = jobPaths.stderr.tailed(100)

  override def preStart() = {
    jobPaths.callRoot.createDirectories()
  }

  lazy val instantiatedScript = {
    def toDockerPath(path: WdlValue): WdlValue = path match {
      case file: WdlFile => WdlFile(jobPaths.toDockerPath(Paths.get(path.valueString)).toAbsolutePath.toString)
      case v => v
    }
    val pathTransformFunction: WdlValue => WdlValue = if (runsOnDocker) toDockerPath else identity

    val localizedInputs = localizeInputs(jobPaths, runsOnDocker, fileSystems, jobDescriptor.inputs)
    call.task.instantiateCommand(localizedInputs, callEngineFunction, pathTransformFunction)
  }

  private def executeScript(script: String): Future[BackendJobExecutionResponse] = {
    logger.info(s"`$script`")
    writeScript(script, if (runsOnDocker) jobPaths.callDockerRoot else jobPaths.callRoot)
    logger.info(s"command: $processArgs")
    process = Option(processArgs.argv.run(ProcessLogger(stdoutWriter writeWithNewline, stderrTailed writeWithNewline)))
    /* This can't run on the context.dispatcher EC because it blocks incoming message processing (in particular abort)
     * TODO use a new separate, backend scoped EC ?
     */
    /*
     * Also note that we only create an asynchronous future now,
     * which guarantees that we won't process another message (say abort) before we have started a process.
     */
    waitAndPostProcess()
  }

  override def execute: Future[BackendJobExecutionResponse] = instantiatedScript match {
    case Success(command) => executeScript(command)
    case Failure(ex) => Future.successful(BackendJobExecutionFailedResponse(jobDescriptor.key, ex))
  }

  override def abortJob: Unit = {
    process foreach { p =>
      p.destroy()
      process = None
    }
  }

  override def recover: Future[BackendJobExecutionResponse] = execute

  /**
    * Writes the script file containing the user's command from the WDL as well
    * as some extra shell code for monitoring jobs
    */
  private def writeScript(instantiatedCommand: String, containerRoot: Path) = {
    jobPaths.script.write(
      s"""#!/bin/sh
          |cd $containerRoot
          |$instantiatedCommand
          |echo $$? > rc
          |""".stripMargin)
  }

  /**
    * --rm automatically deletes the container upon exit
    * -v maps the host workflow executions directory to /root/<workflow id> on the container.
    * -i makes the run interactive, required for the cat and <&0 shenanigans that follow.
    */
  private def buildDockerRunCommand(image: String): String = {
    val dockerDir = jobPaths.callDockerRoot
    s"docker run --rm -v ${jobPaths.callRoot.toAbsolutePath}:$dockerDir -i $image"
  }

  private def waitAndPostProcess(): Future[BackendJobExecutionResponse] = Future {
    val processReturnCode = process map { _.exitValue() } getOrElse (throw new RuntimeException("No process was running"))
    stdoutWriter.writer.flushAndClose()
    stderrTailed.writer.flushAndClose()

    processReturnCode match {
      case ProcessKilledCode => BackendJobExecutionAbortedResponse(jobDescriptor.key) // Special case to check for SIGTERM exit code - implying abort
      case other if other == 0 || runtimeAttributes.dockerImage.isEmpty => postProcessJob()
      case failed =>
        logger.error(s"Non-zero return code: $failed")
        logger.error(s"Standard error was:\n\n${stderrTailed.tailString}\n")
        throw new Exception(s"Unexpected process exit code: $failed")
    }
  }

  private def postProcessJob(): BackendJobExecutionResponse = {
    val stderrFileLength = Try(jobPaths.stderr.size).getOrElse(0L)
    val returnCode = Try(jobPaths.returnCode.contentAsString.stripLineEnd.toInt)
    logger.info(s"Return code: $returnCode")

    if (runtimeAttributes.failOnStderr && stderrFileLength > 0) {
      BackendJobExecutionFailedResponse(jobDescriptor.key, new Throwable(s"Call ${call.fullyQualifiedName}, " +
        s"Workflow ${workflowDescriptor.id}: stderr has length $stderrFileLength"))
    } else {
      lazy val badReturnCodeMessage =
        s"""Call ${call.fullyQualifiedName}, Workflow ${workflowDescriptor.id}: return code was ${returnCode.getOrElse("(none)")}
            |Full command was: $processArgs
            |${stderrTailed.tailString}""".stripMargin

      returnCode match {
        case Success(143) =>
          BackendJobExecutionAbortedResponse(jobDescriptor.key) // Special case to check for SIGTERM exit code - implying abort
        case Success(otherReturnCode) if runtimeAttributes.continueOnReturnCode.continueFor(otherReturnCode) => processSuccess(otherReturnCode)
        case Success(badReturnCode) => BackendJobExecutionFailedResponse(jobDescriptor.key, new Exception(badReturnCodeMessage))
        case Failure(e) => BackendJobExecutionFailedResponse(jobDescriptor.key, new Exception(badReturnCodeMessage, e))
      }
    }
  }

  private def processSuccess(rc: Int) = {
    processOutputs(callEngineFunction, jobPaths) match {
      case Success(outputs) => BackendJobExecutionSucceededResponse(jobDescriptor.key, outputs)
      case Failure(e) =>
        val message = Option(e.getMessage) map { ": " + _ } getOrElse ""
        BackendJobExecutionFailedResponse(jobDescriptor.key, new Throwable("Failed post processing of outputs" + message, e))
    }
  }
}
