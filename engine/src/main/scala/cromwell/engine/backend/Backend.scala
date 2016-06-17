package cromwell.engine.backend

import java.nio.file.{FileSystem, Path}

import akka.actor.ActorSystem
import com.google.api.client.util.ExponentialBackOff
import com.typesafe.config.Config
import cromwell.core.{WorkflowId, WorkflowOptions}
import cromwell.engine.backend.jes.JesBackend
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.backend.runtimeattributes.{ContinueOnReturnCodeFlag, ContinueOnReturnCodeSet, CromwellRuntimeAttributes}
import cromwell.engine.backend.sge.SgeBackend
import cromwell.engine.backend.lsf.LsfBackend
import cromwell.engine.db.DataAccess.ExecutionKeyToJobKey
import cromwell.logging.WorkflowLogger
import cromwell.util.docker.SprayDockerRegistryApiClient
import org.slf4j.LoggerFactory
import wdl4s._
import wdl4s.values._

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Success, Try}

object Backend {
  class StdoutStderrException(message: String) extends RuntimeException(message)

  def from(backendConf: Config, actorSystem: ActorSystem): Backend = Backend.from(backendConf.getString("backend"), actorSystem)

  def from(backendType: BackendType, actorSystem: ActorSystem): Backend = backendType match {
    case BackendType.LOCAL => LocalBackend(actorSystem)
    case BackendType.JES => JesBackend(actorSystem)
    case BackendType.SGE => SgeBackend(actorSystem)
    case BackendType.LSF => LsfBackend(actorSystem)
    case doh => throw new IllegalArgumentException(s"$doh is not a recognized backend")
  }

  def from(name: String, actorSystem: ActorSystem): Backend = {
    val backendType = name.toBackendType
    Backend.from(backendType, actorSystem)
  }

  implicit class BackendyString(val backendType: String) extends AnyVal {
    def toBackendType: BackendType = {
      try {
        BackendType.valueOf(backendType.toUpperCase)
      } catch {
        case e: Exception => throw new IllegalArgumentException(s"$backendType is not a recognized backend")
      }
    }
  }

  private val CallPrefix = "call"
  private val ShardPrefix = "shard"
  private val AttemptPrefix = "attempt"

  def callRootPathWithBaseRoot(jobDescriptor: BackendCallJobDescriptor, baseRoot: String): Path = {
    import io._
    val call = s"$CallPrefix-${jobDescriptor.key.scope.unqualifiedName}"
    val shard = jobDescriptor.key.index map { s => s"$ShardPrefix-$s" } getOrElse ""
    val retry = if (jobDescriptor.key.attempt > 1) s"$AttemptPrefix-${jobDescriptor.key.attempt}" else ""
    val workflowRoot = jobDescriptor.workflowDescriptor.workflowRootPathWithBaseRoot(baseRoot)

    List(call, shard, retry).foldLeft(workflowRoot)((path, dir) => path.resolve(dir)).asDirectory
  }
}

trait BackendJobKey

final case class AttemptedLookupResult(name: String, value: Try[WdlValue]) {
  def toPair = name -> value
}

object AttemptedLookupResult {
  implicit class AugmentedAttemptedLookupSequence(s: Seq[AttemptedLookupResult]) {
    def toLookupMap: Map[String, WdlValue] = s collect {
      case AttemptedLookupResult(name, Success(value)) => (name, value)
    } toMap
  }
}

/**
 * Trait to be implemented by concrete backends.
 */
trait Backend {

  def actorSystem: ActorSystem

  def rootPath(workflowOptions: WorkflowOptions): String

  /**
   * Return a possibly altered copy of inputs reflecting any localization of input file paths that might have
   * been performed for this `Backend` implementation.
   */
  def adjustInputPaths(backendCallJobDescriptor: BackendCallJobDescriptor): CallInputs

  /**
   * Do whatever work is required to initialize the workflow, returning a copy of
   * the coerced inputs present in the `WorkflowDescriptor` with any input `WdlFile`s
   * adjusted for the host workflow execution path.
   */
  def initializeForWorkflow(workflow: WorkflowDescriptor): Future[Unit]

  /**
   * Do whatever cleaning up work is required when a workflow reaches a terminal state.
   */
  def cleanUpForWorkflow(workflow: WorkflowDescriptor)(implicit ec: ExecutionContext): Future[Any] = Future.successful({})

  def engineFunctions(fileSystems: List[FileSystem], workflowContext: WorkflowContext): WorkflowEngineFunctions

  /**
   * Do any work that needs to be done <b>before</b> attempting to restart a workflow.
   */
  def prepareForRestart(restartableWorkflow: WorkflowDescriptor)(implicit ec: ExecutionContext): Future[Unit]

  /**
   * Return CallLogs which contains the stdout/stderr of the particular call
   */
  def stdoutStderr(jobDescriptor: BackendCallJobDescriptor): CallLogs

  def backendType: BackendType

  /**
   * Validate that workflow options contain all required information.
   */
  @throws[IllegalArgumentException]("if a value is missing / incorrect")
  def assertWorkflowOptions(options: WorkflowOptions): Unit = {}

  private[backend] def backendClassString = backendType.toString.toLowerCase.capitalize + "Backend"

  /** Default implementation assumes backends do not support resume, returns an empty Map. */
  def findResumableExecutions(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[ExecutionKeyToJobKey]] = Future.successful(List.empty)

  def workflowLogger(descriptor: WorkflowDescriptor) = WorkflowLogger(
    backendClassString,
    descriptor,
    otherLoggers = Seq(LoggerFactory.getLogger(getClass.getName))
  )

  def jobLogger(jobDescriptor: JobDescriptor[_ <: JobKey]) = WorkflowLogger(
    backendClassString,
    jobDescriptor.workflowDescriptor,
    otherLoggers = Seq(LoggerFactory.getLogger(getClass.getName)),
    callTag = Option(jobDescriptor.key.tag)
  )

  lazy val dockerHashClient = new SprayDockerRegistryApiClient()(actorSystem)

  def pollBackoff: ExponentialBackOff

  def executionInfoKeys: List[String]

  def callRootPathWithBaseRoot(descriptor: BackendCallJobDescriptor, baseRoot: String): Path = Backend.callRootPathWithBaseRoot(descriptor, baseRoot)

  def callRootPath(jobDescriptor: BackendCallJobDescriptor) = {
    callRootPathWithBaseRoot(jobDescriptor, rootPath(jobDescriptor.workflowDescriptor.workflowOptions))
  }

  def callEngineFunctions(descriptor: BackendCallJobDescriptor): CallEngineFunctions

  def instantiateCommand(descriptor: BackendCallJobDescriptor): Try[String]

  def runtimeAttributes(descriptor: BackendCallJobDescriptor): CromwellRuntimeAttributes =
    CromwellRuntimeAttributes(
      descriptor.call.task.runtimeAttributes,
      descriptor,
      Option(descriptor.workflowDescriptor.workflowOptions))

  def poll(jobDescriptor: BackendCallJobDescriptor, previous: ExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle]

  /** Given the specified value for the Docker hash, return the overall hash for this `BackendCall`. */
  private def hashGivenDockerHash(jobDescriptor: BackendCallJobDescriptor)(dockerHash: Option[String]): ExecutionHash = {
    val call = jobDescriptor.call
    val runtimeAttributes = jobDescriptor.callRuntimeAttributes
    val orderedInputs = jobDescriptor.locallyQualifiedInputs.toSeq.sortBy(_._1)
    val orderedOutputs = call.task.outputs.sortWith((l, r) => l.name > r.name)
    val orderedRuntime = Seq(
      ("docker", dockerHash getOrElse ""),
      ("zones", runtimeAttributes.zones.sorted.mkString(",")),
      ("failOnStderr", runtimeAttributes.failOnStderr.toString),
      ("continueOnReturnCode", runtimeAttributes.continueOnReturnCode match {
        case ContinueOnReturnCodeFlag(bool) => bool.toString
        case ContinueOnReturnCodeSet(codes) => codes.toList.sorted.mkString(",")
      }),
      ("cpu", runtimeAttributes.cpu.toString),
      ("preemptible", runtimeAttributes.preemptible.toString),
      ("disks", runtimeAttributes.disks.sortWith((l, r) => l.name > r.name).map(_.toString).mkString(",")),
      ("memoryGB", runtimeAttributes.memoryGB.toString),
      ("bootDiskSizeGb", runtimeAttributes.bootDiskSizeGb.toString)
    )

    val overallHash = Seq(
      backendType.toString,
      call.task.commandTemplateString,
      orderedInputs map { case (k, v) => s"$k=${v.computeHash(jobDescriptor.workflowDescriptor.fileHasher).value}" } mkString "\n",
      orderedRuntime map { case (k, v) => s"$k=$v" } mkString "\n",
      orderedOutputs map { o => s"${o.wdlType.toWdlString} ${o.name} = ${o.requiredExpression.toWdlString}" } mkString "\n"
    ).mkString("\n---\n").md5Sum

    ExecutionHash(overallHash, dockerHash)
  }

  /**
    * Compute a hash that uniquely identifies this call
    */
  def hash(descriptor: BackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[ExecutionHash] = {
    // If a Docker image is defined in the task's runtime attributes, return a `Future[Option[String]]` of the Docker
    // hash string, otherwise return a `Future.successful` of `None`.
    def hashDockerImage(dockerImage: String): Future[Option[String]] =
      if (descriptor.workflowDescriptor.lookupDockerHash)
        dockerHashClient.getDockerHash(dockerImage) map { dh => Option(dh.hashString) }
      else
        Future.successful(Option(dockerImage))

    if (descriptor.workflowDescriptor.writeToCache)
      descriptor.callRuntimeAttributes.docker map hashDockerImage getOrElse Future.successful(None) map hashGivenDockerHash(descriptor)
    else
      Future.successful(ExecutionHash("", None))
  }

  def fileSystems(options: WorkflowOptions): List[FileSystem]

  def buildWorkflowRootPath(rootPath: String, name: String, workflowId: WorkflowId) = s"$rootPath/$name/$workflowId"

  def useCachedCall(cachedCall: BackendCallJobDescriptor, backendCall: BackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[ExecutionHandle]

  def execute(jobDescriptor: BackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[ExecutionHandle]

  def resume(descriptor: BackendCallJobDescriptor, jobKey: BackendJobKey)(implicit ec: ExecutionContext): Future[ExecutionHandle]
}
