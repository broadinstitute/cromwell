package cromwell.backend

import akka.actor.{ActorRef, Props}
import com.typesafe.config.Config
import cromwell.backend.io.WorkflowPathsWithDocker
import cromwell.backend.standard.callcaching.BlacklistCache
import cromwell.core.CallOutputs
import cromwell.core.JobToken.JobTokenType
import cromwell.core.path.Path
import cromwell.core.path.PathFactory.PathBuilders
import cromwell.docker.DockerMirroring
import net.ceedubs.ficus.Ficus._
import wom.expression.{IoFunctionSet, NoIoFunctionSet}
import wom.graph.CommandCallNode

import scala.concurrent.ExecutionContext

trait BackendLifecycleActorFactory extends PlatformSpecific {

  /**
    * Name of the backend.
    *
    * This is the first parameter passed into each factory during creation.
    *
    * @return The configuration-defined name for this instance of the backend.
    */
  def name: String

  def nameForCallCachingPurposes: String =
    configurationDescriptor.backendConfig.getOrElse("name-for-call-caching-purposes", name)

  /**
    * Config values for the backend, and a pointer to the global config.
    *
    * This is the second parameter passed into each factory during creation.
    *
    * @return The backend configuration.
    */
  def configurationDescriptor: BackendConfigurationDescriptor

  /* ****************************** */
  /*     Workflow Initialization    */
  /* ****************************** */

  def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                       ioActor: ActorRef,
                                       calls: Set[CommandCallNode],
                                       serviceRegistryActor: ActorRef,
                                       restarting: Boolean
  ): Option[Props] = None

  /* ****************************** */
  /*          Job Execution         */
  /* ****************************** */

  def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                             initializationData: Option[BackendInitializationData],
                             serviceRegistryActor: ActorRef,
                             ioActor: ActorRef,
                             backendSingletonActor: Option[ActorRef],
                             groupMetricsActor: ActorRef
  ): Props

  lazy val jobExecutionTokenType: JobTokenType = {
    val concurrentJobLimit = configurationDescriptor.backendConfig.as[Option[Int]]("concurrent-job-limit")
    // if defined, use per-backend hog-factor, otherwise use system-level value
    val hogFactor = configurationDescriptor.backendConfig.as[Option[Int]]("hog-factor") match {
      case Some(backendHogFactorValue) => backendHogFactorValue
      case None => configurationDescriptor.globalConfig.getInt("system.hog-safety.hog-factor")
    }
    JobTokenType(name, concurrentJobLimit, hogFactor)
  }

  lazy val jobRestartCheckTokenType: JobTokenType = {
    val concurrentRestartCheckLimit =
      configurationDescriptor.globalConfig.as[Option[Int]]("system.job-restart-check-rate-control.max-jobs")
    // if defined, use per-backend hog-factor, otherwise use system-level value
    val hogFactor = configurationDescriptor.backendConfig.as[Option[Int]]("hog-factor") match {
      case Some(backendHogFactorValue) => backendHogFactorValue
      case None => configurationDescriptor.globalConfig.getInt("system.hog-safety.hog-factor")
    }
    JobTokenType(name, concurrentRestartCheckLimit, hogFactor)
  }

  /* ****************************** */
  /*      Workflow Finalization     */
  /* ****************************** */

  def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                     ioActor: ActorRef,
                                     calls: Set[CommandCallNode],
                                     jobExecutionMap: JobExecutionMap,
                                     workflowOutputs: CallOutputs,
                                     initializationData: Option[BackendInitializationData]
  ): Option[Props] = None

  /* ****************************** */
  /*           Call Caching         */
  /* ****************************** */

  def fileHashingActorProps
    : Option[(BackendJobDescriptor, Option[BackendInitializationData], ActorRef, ActorRef, Option[ActorRef]) => Props] =
    None

  /**
    * Providing this method to generate Props for a cache hit copying actor is optional.
    * To implement it, add a function:
    * def cacheHitCopyingActorInner(jobDescriptor: BackendJobDescriptor,initializationData: Option[BackendInitializationData], serviceRegistryActor: ActorRef, ioActor: ActorRef): Props
    * And then override this method to point to it:
    * override def cacheHitCopyingActorProps = Option(cacheHitCopyingActorInner _)
    *
    * Simples!
    */
  def cacheHitCopyingActorProps: Option[
    (BackendJobDescriptor, Option[BackendInitializationData], ActorRef, ActorRef, Int, Option[BlacklistCache]) => Props
  ] = None

  /* ****************************** */
  /*              Misc.             */
  /* ****************************** */

  def backendSingletonActorProps(serviceRegistryActor: ActorRef): Option[Props] = None

  def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                  jobKey: BackendJobDescriptorKey,
                                  initializationData: Option[BackendInitializationData],
                                  ioActor: ActorRef,
                                  ec: ExecutionContext
  ): IoFunctionSet = NoIoFunctionSet

  def pathBuilders(initializationDataOption: Option[BackendInitializationData]): PathBuilders = List.empty

  def getExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor,
                           backendConfig: Config,
                           initializationData: Option[BackendInitializationData]
  ): Path =
    new WorkflowPathsWithDocker(workflowDescriptor, backendConfig).executionRoot

  def getWorkflowExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor,
                                   backendConfig: Config,
                                   initializationData: Option[BackendInitializationData]
  ): Path =
    new WorkflowPathsWithDocker(workflowDescriptor, backendConfig).workflowRoot

  def runtimeAttributeDefinitions(
    initializationDataOption: Option[BackendInitializationData]
  ): Set[RuntimeAttributeDefinition] = Set.empty

  /**
    * A set of KV store keys that this backend requests that the engine lookup before running each job.
    */
  def requestedKeyValueStoreKeys: Seq[String] = Seq.empty

  /**
    * A set of KV store keys that are requested and looked up on behalf of all backends before running each job.
    */
  def defaultKeyValueStoreKeys: Seq[String] =
    Seq(BackendLifecycleActorFactory.FailedRetryCountKey, BackendLifecycleActorFactory.MemoryMultiplierKey)

  /*
   * Returns credentials that can be used to authenticate to a docker registry server
   * in order to obtain a docker hash.
   */
  def dockerHashCredentials(workflowDescriptor: BackendWorkflowDescriptor,
                            initializationDataOption: Option[BackendInitializationData]
  ): List[Any] = List.empty

  /**
   * Returns a DockerMirror built based on backend configuration
   */
  val dockerMirroring: Option[DockerMirroring] = None

  /**
    * Allows Cromwell to self-identify which cloud it's running on for runtime attribute purposes
    */
  override def platform: Option[Platform] = None
}

object BackendLifecycleActorFactory {
  val FailedRetryCountKey = "FailedRetryCount"
  val MemoryMultiplierKey = "MemoryMultiplier"
}
