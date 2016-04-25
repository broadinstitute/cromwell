package cromwell.backend.impl.local

import java.nio.file.{FileSystems, Path, Paths}

import akka.actor.FSM.NullFunction
import akka.actor.{PoisonPill, Props, LoggingFSM}
import com.typesafe.config.Config
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionAbortedResponse, BackendJobExecutionFailedResponse, BackendJobExecutionResponse, BackendJobExecutionSucceededResponse}
import cromwell.backend.BackendLifecycleActor.{BackendJobExecutionAbortFailedResponse, BackendJobExecutionAbortSucceededResponse, JobAbortResponse}
import cromwell.backend.impl.local.LocalJobExecutionFSM._
import cromwell.core.CallContext
import org.slf4j.LoggerFactory
import wdl4s._
import wdl4s.values.{WdlFile, WdlValue}
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.{Future, Promise}
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Failure, Success, Try}

object LocalJobExecutionFSM {
  sealed trait JobMessages
  case class Run(executePromise: Promise[BackendJobExecutionResponse])
  case class Abort(abortPromise: Promise[JobAbortResponse])
  sealed trait JobState
  case object Idle extends JobState
  case object Running extends JobState
  case object Aborting extends JobState
  protected case class StateData(process: Option[Process])

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

  def props(jobDescriptor: BackendJobDescriptor, backendConfiguration: Config) = Props(new LocalJobExecutionFSM(jobDescriptor, backendConfiguration))
}

class LocalJobExecutionFSM(jobDescriptor: BackendJobDescriptor,
                           backendConfiguration: Config) extends LoggingFSM[JobState, StateData] with SharedFileSystem {
  import LocalJobExecutionActor._
  import better.files._
  import cromwell.core.PathFactory._
  import scala.sys.process._

  startWith(Idle, StateData(None))

  when(Idle) {
    case Event(Run(executePromise), data) => execute(executePromise, data)
    case Event(Abort(abortPromise), _) =>
      abortPromise success BackendJobExecutionAbortFailedResponse(jobDescriptor.key, new RuntimeException("Received Abort message but job was not running."))
      context.stop(self)
      stay()
  }

  when(Running) {
    case Event(Abort(abortPromise), data) => abort(abortPromise, data)
    case Event(Run(_), data) =>
      logger.warn("Received Run message before but job is already running. Ignoring.")
      stay()
  }

  whenUnhandled { NullFunction }

  private def execute(executePromise: Promise[BackendJobExecutionResponse], data: StateData) = {
    jobPaths.callRoot.createDirectories()

    instantiatedScript match {
      case Success(command) =>
        val process = spawnProcess(command)
        val futureResult = waitAndPostProcess(process)
        futureResult onComplete {
          case _ => self ! PoisonPill
        }
        executePromise completeWith futureResult
        goto(Running) using data.copy(process = Option(process))
      case Failure(ex) =>
        executePromise success BackendJobExecutionFailedResponse(jobDescriptor.key, ex)
        context.stop(self)
        stay()
    }
  }

  private def abort(abortPromise: Promise[JobAbortResponse], data: StateData) = {
    data.process map { p =>
      p.destroy()
      abortPromise success BackendJobExecutionAbortSucceededResponse(jobDescriptor.key)
    } getOrElse {
      abortPromise success BackendJobExecutionAbortFailedResponse(jobDescriptor.key, new RuntimeException(s"Tried to abort ${jobDescriptor.key} before executing it."))
    }
    context.stop(self)
    stay()
  }

  val workflowDescriptor = jobDescriptor.descriptor
  val jobPaths = new JobPaths(workflowDescriptor, backendConfiguration, jobDescriptor.key)
  val fileSystemsConfig = backendConfiguration.getConfig("filesystems")
  override val sharedFsConfig = fileSystemsConfig.getConfig("local")

  val call = jobDescriptor.key.call
  val callEngineFunction = {
    val callContext = new CallContext(
      jobPaths.callRoot.toString,
      jobPaths.stdout.toAbsolutePath.toString,
      jobPaths.stderr.toAbsolutePath.toString
    )

    new LocalCallEngineFunctions(fileSystems, callContext)
  }

  val lookup = {
    val declarations = workflowDescriptor.workflowNamespace.workflow.declarations ++ call.task.declarations
    val unqualifiedWorkflowInputs = workflowDescriptor.inputs map {
      case (fqn, v) => splitFqn(fqn)._2 -> v
    }
    val knownInputs = unqualifiedWorkflowInputs ++ jobDescriptor.symbolMap
    WdlExpression.standardLookupFunction(knownInputs, declarations, callEngineFunction)
  }

  private def evaluate(wdlExpression: WdlExpression) = wdlExpression.evaluate(lookup, callEngineFunction)

  val runtimeAttributes = {
    val evaluateAttrs = call.task.runtimeAttributes.attrs mapValues evaluate
    // Fail the call if runtime attributes can't be evaluated
    val evaluatedAttributes = TryUtils.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
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

  lazy val instantiatedScript = {
    def toDockerPath(path: WdlValue): WdlValue = path match {
      case file: WdlFile => WdlFile(jobPaths.toDockerPath(Paths.get(path.valueString)).toAbsolutePath.toString)
      case v => v
    }
    val pathTransformFunction: WdlValue => WdlValue = if (runsOnDocker) toDockerPath else identity

    // Inputs coming from the workflow inputs (json input mapping)
    val workflowInputEntries = workflowDescriptor.inputs collect {
      case (fqn, value) if splitFqn(fqn)._1 == call.fullyQualifiedName => splitFqn(fqn)._2 -> value
    }

    // Inputs coming from the "input" keyword in the workflow declaration. These need to be evaluated because they're expressions
    val evaluatedInputMappings = call.inputMappings mapValues evaluate

    TryUtils.sequenceMap(evaluatedInputMappings, "Job Input evaluation") flatMap { inputs =>
      val localizedInputs = localizeInputs(jobPaths, runsOnDocker, fileSystems, inputs ++ workflowInputEntries)
      call.task.instantiateCommand(localizedInputs, callEngineFunction, pathTransformFunction)
    }
  }

  private def spawnProcess(script: String): Process = {
    logger.info(s"`$script`")
    writeScript(script, if (runsOnDocker) jobPaths.callDockerRoot else jobPaths.callRoot)
    logger.info(s"command: $processArgs")
    processArgs.argv.run(ProcessLogger(stdoutWriter writeWithNewline, stderrTailed writeWithNewline))
  }

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

  private def waitAndPostProcess(process: Process): Future[BackendJobExecutionResponse] = Future {
    val processReturnCode = process.exitValue()
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
    processOutputs(jobDescriptor, workflowDescriptor.id, lookup, callEngineFunction, jobPaths) match {
      case Success(outputs) => BackendJobExecutionSucceededResponse(jobDescriptor.key, outputs)
      case Failure(e) =>
        val message = Option(e.getMessage) map { ": " + _ } getOrElse ""
        BackendJobExecutionFailedResponse(jobDescriptor.key, new Throwable("Failed post processing of outputs" + message, e))
    }
  }

}
