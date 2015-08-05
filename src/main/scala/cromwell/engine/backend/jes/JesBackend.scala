package cromwell.engine.backend.jes

import java.math.BigInteger
import java.nio.file.{Path, Paths}

import com.google.api.services.genomics.model.Parameter
import com.typesafe.config.ConfigFactory
import com.typesafe.scalalogging.LazyLogging
import cromwell.binding.WdlExpression._
import cromwell.binding._
import cromwell.binding.types.WdlFileType
import cromwell.binding.values._
import cromwell.engine.WorkflowId
import cromwell.engine.backend.{StdoutStderr, Backend}
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.db.DataAccess
import cromwell.parser.BackendType
import cromwell.util.TryUtil
import cromwell.util.google.GoogleCloudStoragePath

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

object JesBackend {
  private lazy val JesConf = ConfigFactory.load.getConfig("backend").getConfig("jes")
  lazy val GoogleProject = JesConf.getString("project")
  lazy val GoogleApplicationName = JesConf.getString("applicationName")
  lazy val CromwellExecutionBucket = JesConf.getString("baseExecutionBucket")

  lazy val JesConnection = JesInterface(GoogleApplicationName)

  /*
    FIXME: At least for now the only files that can be used are stdout/stderr. However this leads to a problem
    where stdout.txt is input and output. Redirect stdout/stderr to a different name, but it'll be localized back
    in GCS as stdout/stderr. Yes, it's hacky.
   */
  val LocalStdoutParamName = "job_stdout"
  val LocalStderrParamName = "job_stderr"

  val LocalStdoutValue = "job.stdout.txt"
  val LocalStderrValue = "job.stderr.txt"

  private def callGcsPath(workflowId: String, workflowName: String, callName: String): String =
    s"$CromwellExecutionBucket/$workflowName/$workflowId/call-$callName"

  // Decoration around WorkflowDescriptor to generate bucket names and the like
  implicit class JesWorkflowDescriptor(val descriptor: WorkflowDescriptor) extends AnyVal {
    def callDir(call: Call) = callGcsPath(descriptor.id.toString, descriptor.name, call.name)
  }

  def stderrJesOutput(callGcsPath: String): JesOutput = JesOutput(LocalStderrParamName, s"$callGcsPath/$LocalStderrValue", Paths.get(LocalStderrValue))
  def stdoutJesOutput(callGcsPath: String): JesOutput = JesOutput(LocalStdoutParamName, s"$callGcsPath/$LocalStdoutValue", Paths.get(LocalStdoutValue))

  // For now we want to always redirect stdout and stderr. This could be problematic if that's what the WDL calls stuff, but oh well
  def standardParameters(callGcsPath: String): Seq[JesParameter] = Seq(
    stdoutJesOutput(callGcsPath),
    stderrJesOutput(callGcsPath)
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

  /**
   * Takes a path in GCS and comes up with a local path which is unique for the given GCS path
   * @param gcsPath The input path
   * @return A path which is unique per input path
   */
  def localFilePathFromCloudStoragePath(gcsPath: GoogleCloudStoragePath): Path = {
    Paths.get("/cromwell_root/" + gcsPath.bucket + "/" + gcsPath.objectName)
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
  override def initializeForWorkflow(workflow: WorkflowDescriptor): HostInputs = workflow.actualInputs

  def taskOutputToJesOutput(taskOutput: TaskOutput, callGcsPath: String, scopedLookupFunction: ScopedLookupFunction, engineFunctions: JesEngineFunctions): Try[Option[JesOutput]] = {
    taskOutput.wdlType match {
      case WdlFileType =>
        taskOutput.expression.evaluate(scopedLookupFunction, engineFunctions, interpolateStrings = true) match {
          case Success(v) => Success(Option(JesOutput(taskOutput.name, s"$callGcsPath/${taskOutput.name}", Paths.get(v.valueString))))
          case Failure(e) => Failure(new IllegalArgumentException(s"JES requires File outputs to be determined prior to running, but ${taskOutput.name} can not."))
        }
      case _ => Success(None)
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

  override def stdoutStderr(workflowId: WorkflowId, workflowName: String, callName: String): StdoutStderr = {
    val base = callGcsPath(workflowId.toString, workflowName, callName)
    StdoutStderr(
      stdout = WdlFile(s"$base/$LocalStdoutValue"),
      stderr = WdlFile(s"$base/$LocalStderrValue")
    )
  }

  override def executeCommand(instantiatedCommandLine: String,
                              workflowDescriptor: WorkflowDescriptor,
                              call: Call,
                              backendInputs: CallInputs,
                              scopedLookupFunction: ScopedLookupFunction): Try[Map[String, WdlValue]] = {
    val callGcsPath = s"${workflowDescriptor.callDir(call)}"

    val engineFunctions = new JesEngineFunctions(GoogleCloudStoragePath(callGcsPath), JesConnection)

    // FIXME: Not particularly robust at the moment.
    val jesInputs: Seq[JesParameter] = backendInputs.collect({
      case (k, v) if v.isInstanceOf[WdlFile] => JesInput(k, scopedLookupFunction(k).valueString, Paths.get(v.valueString))
    }).toSeq

    // Call preevaluateExpressionForFilenames for each of the task output expressions, and flatten the lists into a single Try[Seq[WdlFile]]
    val filesWithinExpressions = Try(
      call.task.outputs.map {
        taskOutput => taskOutput.expression.preevaluateExpressionForFilenames(scopedLookupFunction, engineFunctions)
      }.flatMap({_.get})
    )

    val jesOutputs: Try[Seq[JesParameter]] = filesWithinExpressions match {
      case Failure(error) => Failure(error)
      case Success(fileSeq) =>
        val anonOutputs: Seq[JesParameter] = fileSeq map {x => anonymousTaskOutput(x.value, engineFunctions)}
        // FIXME: If localizeTaskOutputs gives a Failure, need to Fail entire function - don't use .get
        val namedOutputs: Seq[JesParameter] = taskOutputsToJesOutputs(call.task.outputs, callGcsPath, scopedLookupFunction, engineFunctions).get
        Success(anonOutputs ++ namedOutputs)
    }

    // FIXME: Ignore all the errors!
    val unsafeJesOutputs: Seq[JesParameter] = jesOutputs.get

    val jesParameters = standardParameters(callGcsPath) ++ jesInputs ++ unsafeJesOutputs

    // Not possible to currently get stdout/stderr so redirect everything and hope the WDL isn't doing that too
    val redirectedCommand = s"$instantiatedCommandLine > $LocalStdoutValue 2> $LocalStderrValue"

    val status = Pipeline(redirectedCommand, workflowDescriptor, call, jesParameters, GoogleProject, JesConnection).run.waitUntilComplete()

    def taskOutputToRawValue(taskOutput: TaskOutput): Try[WdlValue] = {
      val jesOutputFile = unsafeJesOutputs find {_.name == taskOutput.name} map {j => WdlFile(j.gcs)}
      jesOutputFile.map(Success(_)).getOrElse(taskOutput.expression.evaluate(scopedLookupFunction, engineFunctions, interpolateStrings = true))
    }

    val outputMappings = call.task.outputs.map { taskOutput =>
      val rawValue = taskOutputToRawValue(taskOutput)
      logger.debug(s"JesBackend setting ${taskOutput.name} to $rawValue")
      taskOutput.name -> rawValue
    }.toMap

    val warnAboutStderrLength: BigInteger =
      if (!call.failOnStderr) BigInteger.valueOf(0)
      else JesConnection.storage.objectSize(GoogleCloudStoragePath(stderrJesOutput(callGcsPath).gcs))

    if (warnAboutStderrLength.intValue > 0) {
      Failure(new Throwable(s"Workflow ${workflowDescriptor.id}: stderr has length $warnAboutStderrLength for command: $instantiatedCommandLine"))
    } else {
      status match {
        //
        case Run.Success(created, started, finished) =>
          unwrapOutputValues(outputMappings, workflowDescriptor)
        case Run.Failed(created, started, finished, errorCode, errorMessage) =>
          Failure(new Throwable(s"Workflow ${workflowDescriptor.id}: errorCode $errorCode for command: $instantiatedCommandLine. Message: $errorMessage"))
      }
    }
  }

  def unwrapOutputValues(outputMappings: Map[String, Try[WdlValue]], workflowDescriptor: WorkflowDescriptor): Try[Map[String, WdlValue]] = {
    val taskOutputEvaluationFailures = outputMappings.filter { _._2.isFailure }
    if (taskOutputEvaluationFailures.isEmpty) {
      Success(outputMappings.collect { case (name, Success(wdlValue)) => name -> wdlValue })
    } else {
      val message = taskOutputEvaluationFailures.collect { case (name, Failure(e)) => s"$name: $e" }.mkString("\n")
      Failure(new Throwable(s"Workflow ${workflowDescriptor.id}: $message"))
    }
  }

  override def handleCallRestarts(restartableWorkflows: Seq[RestartableWorkflow], dataAccess: DataAccess)(implicit ec: ExecutionContext): Future[Any] = Future("FIXME")

  override def backendType = BackendType.JES
}
