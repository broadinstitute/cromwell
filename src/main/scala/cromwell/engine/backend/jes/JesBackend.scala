package cromwell.engine.backend.jes

import java.io.File
import java.math.BigInteger
import java.net.URL
import java.nio.file.{Path, Paths}

import com.google.api.services.genomics.model.Parameter
import com.typesafe.config.ConfigFactory
import com.typesafe.scalalogging.LazyLogging
import cromwell.binding.WdlExpression._
import cromwell.binding._
import cromwell.binding.types.WdlFileType
import cromwell.binding.values._
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortFunction, AbortRegistrationFunction}
import cromwell.engine.backend.{BackendCall, TaskAbortedException, Backend, StdoutStderr}
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.db.DataAccess
import cromwell.engine.WorkflowId
import cromwell.parser.BackendType
import cromwell.util.TryUtil
import cromwell.util.google.GoogleCloudStoragePath

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}
import scala.language.postfixOps

object JesBackend {
  private lazy val JesConf = ConfigFactory.load.getConfig("backend").getConfig("jes")
  lazy val GoogleProject = JesConf.getString("project")
  lazy val GoogleApplicationName = JesConf.getString("applicationName")
  lazy val CromwellExecutionBucket = JesConf.getString("baseExecutionBucket")
  lazy val EndpointUrl = new URL(JesConf.getString("endpointUrl"))
  lazy val JesConnection = JesInterface(GoogleApplicationName, EndpointUrl)

  /*
    FIXME: At least for now the only files that can be used are stdout/stderr. However this leads to a problem
    where stdout.txt is input and output. Redirect stdout/stderr to a different name, but it'll be localized back
    in GCS as stdout/stderr. Yes, it's hacky.
   */
  val LocalStdoutParamName = "job_stdout"
  val LocalStderrParamName = "job_stderr"

  val LocalStdoutValue = "job.stdout.txt"
  val LocalStderrValue = "job.stderr.txt"

  val JesCromwellRoot = "/cromwell_root"

  def callGcsPath(workflowId: String, workflowName: String, callName: String, index: ExecutionIndex): String = {
    val shardPath = index map { i => s"/shard-$i"} getOrElse ""
    s"$CromwellExecutionBucket/$workflowName/$workflowId/call-$callName$shardPath"
  }

  // Decoration around WorkflowDescriptor to generate bucket names and the like
  implicit class JesWorkflowDescriptor(val descriptor: WorkflowDescriptor) extends AnyVal {
    def callDir(key: CallKey) = callGcsPath(descriptor.id.toString, descriptor.name, key.scope.name, key.index)
  }

  def stderrJesOutput(callGcsPath: String): JesOutput = JesOutput(LocalStderrParamName, s"$callGcsPath/$LocalStderrValue", Paths.get(LocalStderrValue))
  def stdoutJesOutput(callGcsPath: String): JesOutput = JesOutput(LocalStdoutParamName, s"$callGcsPath/$LocalStdoutValue", Paths.get(LocalStdoutValue))
  def localizationDiskInput(): JesInput = JesInput("working_disk", "disk://local-disk", new File(JesCromwellRoot).toPath)

  // For now we want to always redirect stdout and stderr. This could be problematic if that's what the WDL calls stuff, but oh well
  def standardParameters(callGcsPath: String): Seq[JesParameter] = Seq(
    stdoutJesOutput(callGcsPath),
    stderrJesOutput(callGcsPath),
    localizationDiskInput()
  )

  sealed trait JesParameter {
    def name: String
    def gcs: String
    def local: Path

    final val isInput = this.isInstanceOf[JesInput]
    final val isOutput = !isInput
    final val toGoogleParameter = new Parameter().setName(name).setValue(local.toString).setType("REFERENCE")
  }

  final case class JesInput(name: String, gcs: String, local: Path) extends JesParameter
  final case class JesOutput(name: String, gcs: String, local: Path) extends JesParameter
}

class JesBackend extends Backend with LazyLogging {
  type BackendCall = JesBackendCall

  /**
   * Takes a path in GCS and comes up with a local path which is unique for the given GCS path
   * @param gcsPath The input path
   * @return A path which is unique per input path
   */
  def localFilePathFromCloudStoragePath(gcsPath: GoogleCloudStoragePath): Path = {
    Paths.get(JesCromwellRoot + "/" + gcsPath.bucket + "/" + gcsPath.objectName)
  }

  /**
   * Takes a single FQN to WdlValue pair and maps google cloud storage (GCS) paths into an appropriate local file path.
   * If the input is not a WdlFile, or the WdlFile is not a GCS path, the mapping is a noop.
   *
   * @param fqn the FQN of the input variable
   * @param wdlValue the value of the input
   * @return a new FQN to WdlValue pair, with WdlFile paths modified if appropriate.
   */
  def mapInputValue(fqn: String, wdlValue: WdlValue): (String, WdlValue) = {
    wdlValue match {
      case WdlFile(path) =>
        GoogleCloudStoragePath.parse(path) match {
          case Success(gcsPath) => (fqn, WdlFile(localFilePathFromCloudStoragePath(gcsPath).toString))
          case Failure(e) => (fqn, wdlValue)
        }
      case array:WdlArray => fqn -> array.map(v => mapInputValue("bogus", v)._2) // TODO: this is ugly.
      case _ => (fqn, wdlValue)
    }
  }

  override def adjustInputPaths(call: Call, inputs: CallInputs): CallInputs = inputs map {case (k,v) => mapInputValue(k,v)}
  override def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs = outputs

  // No need to copy GCS inputs for the workflow we should be able to directly reference them
  override def initializeForWorkflow(workflow: WorkflowDescriptor): Try[HostInputs] = Success(workflow.actualInputs)

  def taskOutputToJesOutput(taskOutput: TaskOutput, callGcsPath: String, scopedLookupFunction: ScopedLookupFunction, engineFunctions: JesEngineFunctions): Try[Option[JesOutput]] = {
    // If the output isn't a file then it's not a file to be retrieved by the JES run! (but NB: see anonymous task output below)
    // Currently, no function calls can be run before the JES call completes, so we can't retrieve files with names based on function calls.
    if (taskOutput.wdlType != WdlFileType || taskOutput.expression.containsFunctionCall) {
      Success(None)
    } else {
      val evaluatedExpression = taskOutput.expression.evaluate(scopedLookupFunction, engineFunctions, interpolateStrings = true)
      evaluatedExpression match {
        case Success(v) =>
          val jesOutput = JesOutput(taskOutput.name, s"$callGcsPath/${taskOutput.name}", Paths.get(v.valueString))
          Success(Option(jesOutput))
        case Failure(e) => Failure(new IllegalArgumentException(s"JES requires File outputs to be determined prior to running, but ${taskOutput.name} can not."))
      }
    }
  }

  def taskOutputsToJesOutputs(taskOutputs: Seq[TaskOutput], callGcsPath: String, scopedLookupFunction: ScopedLookupFunction, engineFunctions: JesEngineFunctions): Try[Seq[JesOutput]] = {
    val localizedOutputs = taskOutputs map {taskOutputToJesOutput(_, callGcsPath, scopedLookupFunction, engineFunctions)}
    val localizationFailures = localizedOutputs filter {_.isFailure}
    if (localizationFailures.isEmpty) {
      Success(localizedOutputs.collect({case Success(x) => x}).flatten)
    } else {
      val failureMessages = TryUtil.stringifyFailures(localizationFailures)
      Failure(new IllegalArgumentException(failureMessages.mkString("\n")))
    }
  }

  def anonymousTaskOutput(value: String, engineFunctions: JesEngineFunctions): JesOutput = {
    JesOutput(value, engineFunctions.gcsPathFromAnyString(value).toString, Paths.get(value))
  }

  override def stdoutStderr(workflowId: WorkflowId, workflowName: String, callName: String, index: ExecutionIndex): StdoutStderr = {
    val base = callGcsPath(workflowId.toString, workflowName, callName, index)
    StdoutStderr(
      stdout = WdlFile(s"$base/$LocalStdoutValue"),
      stderr = WdlFile(s"$base/$LocalStderrValue")
    )
  }

  override def bindCall(workflowDescriptor: WorkflowDescriptor,
                        key: CallKey,
                        locallyQualifiedInputs: CallInputs,
                        abortRegistrationFunction: AbortRegistrationFunction): BackendCall = {
    JesBackendCall(this, workflowDescriptor, key, locallyQualifiedInputs, abortRegistrationFunction)
  }

  def execute(backendCall: BackendCall): Try[CallOutputs] = {
    val tag = makeTag(backendCall)
    val cmdInput = JesInput("exec", backendCall.gcsExecPath.toString, Paths.get("exec.sh"))
    logger.info(s"$tag Call GCS path: ${backendCall.callGcsPath}")

    // FIXME: Not particularly robust at the moment.
    val jesInputs: Seq[JesParameter] = adjustInputPaths(backendCall.call, backendCall.locallyQualifiedInputs).collect({
      case (k, v) if v.isInstanceOf[WdlFile] =>
        logger.info(s"$tag: $k -> @${v.valueString}@")
        JesInput(k, backendCall.lookupFunction(k).valueString, Paths.get(v.valueString))
    }).toSeq :+ cmdInput

    val filesWithinExpressions = Try(
      backendCall.call.task.outputs map {
        taskOutput => taskOutput.expression.preevaluateExpressionForFilenames(backendCall.lookupFunction, backendCall.engineFunctions)
      } flatMap { _.get }
    )

    val jesOutputs: Try[Seq[JesParameter]] = filesWithinExpressions match {
      case Failure(error) => Failure(error)
      case Success(fileSeq) =>
        val anonOutputs: Seq[JesParameter] = fileSeq map {x => anonymousTaskOutput(x.value, backendCall.engineFunctions)}
        // FIXME: If localizeTaskOutputs gives a Failure, need to Fail entire function - don't use .get
        val namedOutputs: Seq[JesParameter] = taskOutputsToJesOutputs(
          backendCall.call.task.outputs,
          backendCall.callGcsPath,
          backendCall.lookupFunction,
          backendCall.engineFunctions
        ).get
        Success(anonOutputs ++ namedOutputs)
    }

    backendCall.instantiateCommand match {
      case Success(command) => runWithJes(backendCall, command, jesInputs, jesOutputs.get)
      case Failure(ex) => Failure(ex)
    }
  }

  private def runWithJes(backendCall: BackendCall, command: String, jesInputs: Seq[JesParameter], jesOutputs: Seq[JesParameter]): Try[CallOutputs] = {
    val tag = makeTag(backendCall)
    // FIXME: Ignore all the errors!
    val unsafeJesOutputs: Seq[JesParameter] = jesOutputs
    val jesParameters = standardParameters(backendCall.callGcsPath) ++ jesInputs ++ unsafeJesOutputs
    logger.info(s"$tag `$command`")
    JesConnection.storage.uploadObject(backendCall.gcsExecPath, command)

    val run = Pipeline(s"/bin/bash exec.sh > $LocalStdoutValue 2> $LocalStderrValue", backendCall.workflowDescriptor, backendCall.key, jesParameters, GoogleProject, JesConnection).run
    // Wait until the job starts (or completes/fails) before registering the abort to avoid awkward cancel-during-initialization behavior.
    run.waitUntilRunningOrComplete()
    backendCall.callAbortRegistrationFunction.register(AbortFunction(() => run.abort()))

    try {
      val status = run.waitUntilComplete(None)

      def taskOutputToRawValue(taskOutput: TaskOutput): Try[WdlValue] = {
        val jesOutputFile = unsafeJesOutputs find { _.name == taskOutput.name } map { j => WdlFile(j.gcs) }
        jesOutputFile.map(Success(_)).getOrElse(taskOutput.expression.evaluate(backendCall.lookupFunction, backendCall.engineFunctions, interpolateStrings = true))
      }

      def processOutput(taskOutput: TaskOutput, processingFunction: TaskOutput => Try[WdlValue]) = {
        val rawValue = processingFunction(taskOutput)
        logger.debug(s"JesBackend setting ${taskOutput.name} to $rawValue")
        taskOutput.name -> rawValue
      }

      val outputMappings = backendCall.call.task.outputs map {taskOutput => processOutput(taskOutput, taskOutputToRawValue)} toMap

      lazy val stderrLength: BigInteger = JesConnection.storage.objectSize(GoogleCloudStoragePath(stderrJesOutput(backendCall.callGcsPath).gcs))

      if (backendCall.call.failOnStderr && stderrLength.intValue > 0) {
        Failure(new Throwable(s"Workflow ${backendCall.workflowDescriptor.id}: stderr has length $stderrLength for command: $command"))
      } else status match {
        case Run.Success =>
          unwrapOutputValues(outputMappings, backendCall.workflowDescriptor)
        case Run.Failed(errorCode, errorMessage) =>
          val throwable = if (errorMessage contains "Operation canceled at") {
            new TaskAbortedException()
          } else {
            new Throwable(s"Workflow ${backendCall.workflowDescriptor.id}: errorCode $errorCode for command: $command. Message: $errorMessage")
          }
          Failure(throwable)
      }
    }
    catch {
      case e: Exception => Failure(e)
    }
  }

  def unwrapOutputValues(outputMappings: Map[String, Try[WdlValue]], workflowDescriptor: WorkflowDescriptor): Try[Map[String, WdlValue]] = {
    val taskOutputEvaluationFailures = outputMappings filter { _._2.isFailure }
    if (taskOutputEvaluationFailures.isEmpty) {
      Success(outputMappings collect { case (name, Success(wdlValue)) => name -> wdlValue })
    } else {
      val message = taskOutputEvaluationFailures collect { case (name, Failure(e)) => s"$name: $e" } mkString "\n"
      Failure(new Throwable(s"Workflow ${workflowDescriptor.id}: $message"))
    }
  }

  override def handleCallRestarts(restartableWorkflows: Seq[RestartableWorkflow], dataAccess: DataAccess)(implicit ec: ExecutionContext): Future[Any] = Future("FIXME")

  override def backendType = BackendType.JES
}
