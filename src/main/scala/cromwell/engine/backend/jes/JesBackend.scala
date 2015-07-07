package cromwell.engine.backend.jes

import java.nio.file.{Path, Paths}

import com.google.api.services.genomics.model.{Parameter, ServiceAccount}
import com.typesafe.config.ConfigFactory
import com.typesafe.scalalogging.LazyLogging
import cromwell.binding.WdlExpression._
import cromwell.binding._
import cromwell.binding.values._
import cromwell.engine.backend.Backend
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.db.DataAccess
import cromwell.util.TryUtil
import cromwell.util.google.{GoogleScopes, GoogleCloudStoragePath}
import scala.collection.JavaConverters._
import scala.concurrent.{Future, ExecutionContext}
import scala.util.{Failure, Success, Try}
import JesBackend._

object JesBackend {
  private lazy val JesConf = ConfigFactory.load.getConfig("backend").getConfig("jes")
  lazy val GoogleProject = JesConf.getString("project")
  lazy val GoogleApplicationName = JesConf.getString("applicationName")

  lazy val JesConnection = JesInterface(GoogleApplicationName)

  /*
    FIXME: At least for now the only files that can be used are stdout/stderr. However this leads to a problem
    where stdout.txt is input and output :) Redirect stdout/stderr to a different name, but it'll be localized back
    in GCS as stdout/stderr. Yes, it's hacky.
   */
  val LocalStdout = "job.stdout.txt"
  val LocalStderr = "job.stderr.txt"

  val CromwellExecutionBucket = s"gs://cromwell-dev/cromwell-executions"

  def gsInputToLocal(gsPath: String): Path = {
    // FIXME: What if it doesn't have a gs?
    Paths.get(gsPath.replaceFirst("^gs:/", "/tmp"))
  }

  // Decoration around WorkflowDescriptor to generate bucket names and the like
  implicit class JesWorkflowDescriptor(val descriptor: WorkflowDescriptor) extends AnyVal {
    def bucket = s"$CromwellExecutionBucket/${descriptor.name}/${descriptor.id.toString}"
    def callDir(call: Call) = s"$bucket/call-${call.name}"
  }

  // For now we want to always redirect stdout and stderr. This could be problematic if that's what the WDL calls stuff, but oh well
  def standardParameters(callGcsPath: String): Seq[JesParameter] = Seq(
    JesOutput(LocalStderr, s"$callGcsPath/stderr.txt", Paths.get("stderr.txt")),
    JesOutput(LocalStdout, s"$callGcsPath/stdout.txt", Paths.get("stdout.txt"))
  )

  sealed trait JesParameter { // FIXME: Perhaps not the best name
    def name: String
    def gcs: String
    def local: Path

    final def isInput = this.isInstanceOf[JesInput]
    final def isOutput = !isInput
    // FIXME: Perhaps not the best name
    final def toGoogleParamter = new Parameter().setName(name).setValue(local.toString).setType("REFERENCE")
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
    Paths.get("/some_unlikely_folder/" + gcsPath.bucket + "/" + gcsPath.objectName)
  }

  /**
   * Takes a single input mapping from FQN to WdlValue and maps google cloud storage (GCS) paths into an appropriate local file path.
   * If the input is not a WdlFile, or the WdlFile is not a GCS path, the mapping is a noop.
   *
   * @param fqn the FQN of the input variable
   * @param wdlValue the value of the input
   * @return a new FQN to WdlValue pair, with WdlFile paths modified if appropriate.
   */
  def mapInputValue(fqn: String, wdlValue: WdlValue): (String, WdlValue) = {
    wdlValue match {
      case WdlFile(path) => {
        GoogleCloudStoragePath.parse(path) match {
          case Success(gcsPath) => (fqn, WdlFile(localFilePathFromCloudStoragePath(gcsPath).toString))
          case Failure(e) => (fqn, wdlValue)
        }
      }
      case _ => (fqn, wdlValue)
    }
  }

  override def adjustInputPaths(call: Call, inputs: CallInputs): CallInputs = {
    inputs map { case (k,v) => mapInputValue(k,v) }
  }

  override def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs = outputs
  
  // No need to copy GCS inputs for the workflow we should be able to direclty reference them
  override def initializeForWorkflow(workflow: WorkflowDescriptor): HostInputs = workflow.actualInputs

  // FIXME: signature is weird
  def localizeTaskOutput(taskOutput: TaskOutput, callGcsPath: String, scopedLookupFunction: ScopedLookupFunction, engineFunctions: JesEngineFunctions): Try[Option[JesOutput]] = {
    taskOutput.wdlType match {
      case f: WdlFile =>
        taskOutput.expression.evaluate(scopedLookupFunction, engineFunctions) match {
          case Success(v) => Success(Option(JesOutput(taskOutput.name, s"$callGcsPath/${taskOutput.name}", Paths.get(v.valueString))))
          case Failure(e) => Failure(new IllegalArgumentException(s"JES requires File outputs to be determined prior to running, but ${taskOutput.name} can not."))
        }
      case _ => Success(None)
    }
  }

  // FIXME: Signature is weird
  def localizeTaskOutputs(taskOutputs: Seq[TaskOutput], callGcsPath: String, scopedLookupFunction: ScopedLookupFunction, engineFunctions: JesEngineFunctions): Try[Seq[JesOutput]] = {
    val localizedOutputs = taskOutputs map {localizeTaskOutput(_, callGcsPath, scopedLookupFunction, engineFunctions)}
    val localizationFailures = localizedOutputs filter {_.isFailure}
    if (localizationFailures.isEmpty) {
      Success(localizedOutputs.collect({case Success(x) => x}).flatten)
    } else {
      val failureMessages = TryUtil.stringifyFailures(localizationFailures)
      Failure(new IllegalArgumentException(failureMessages.mkString("\n")))
    }
  }
}

class JesBackend extends Backend with LazyLogging {
  // FIXME: As of right now this isn't needed, but it probably will be as things get less hacky
  override def adjustInputPaths(call: Call, inputs: CallInputs): CallInputs = inputs

  // FIXME: Pass through for now, but perhaps not as time goes on
  override def initializeForWorkflow(workflow: WorkflowDescriptor): HostInputs = {
    // FIXME: No need to copy GCS inputs for the workflow we should be able to direclty reference them
    workflow.actualInputs
  }

  override def executeCommand(commandLine: String,
                              workflowDescriptor: WorkflowDescriptor,
                              call: Call,
                              backendInputs: CallInputs,
                              scopedLookupFunction: ScopedLookupFunction):Try[Map[String, WdlValue]] = ???
//  {
//    val gcsPath = s"${workflowDescriptor.callDir(call)}"
//
//    // FIXME: These are only used for that TemporaryEngineFunctions at the moment
//    val stdout = s"$gcsPath/stdout.txt"
//    val stderr = s"$gcsPath/stderr.txt"
//
//    // FIXME: Not possible to currently get stdout/stderr so redirect everything and hope the WDL isn't doing that too
//    val redirectedCommand = s"$commandLine > $LocalStdout  2> $LocalStderr"
//
//    Pipeline(redirectedCommand, workflowDescriptor, call, backendInputs, GoogleProject, GenomicsService).run.waitUntilComplete()
//
//    // FIXME: This is *only* going to work w/ my hacked stdout/stderr above, see LocalBackend for more info
//    call.task.outputs.map { taskOutput =>
//      val rawValue = taskOutput.expression.evaluate(
//        scopedLookupFunction,
//        new TemporaryEngineFunctions(gcsPath, stdout, stderr)
//      )
//      println(s"JesBackend setting ${taskOutput.name} to $rawValue")
//      taskOutput.name -> rawValue
//    }.toMap
//  }

  override def handleCallRestarts(restartableWorkflows: Seq[RestartableWorkflow], dataAccess: DataAccess)(implicit ec: ExecutionContext): Future[Any] = ???
}
