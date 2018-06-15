package cromwell.server

import akka.actor.SupervisorStrategy.{Escalate, Restart}
import akka.actor.{Actor, ActorInitializationException, ActorLogging, ActorRef, OneForOneStrategy}
import akka.event.Logging
import akka.http.scaladsl.Http
import akka.pattern.GracefulStopSupport
import akka.routing.RoundRobinPool
import akka.stream.ActorMaterializer
import com.typesafe.config.ConfigFactory
import cromwell.core._
import cromwell.core.actor.StreamActorHelper.ActorRestartException
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.core.io.Throttle
import cromwell.docker.DockerHashActor
import cromwell.docker.DockerHashActor.DockerHashContext
import cromwell.docker.local.DockerCliFlow
import cromwell.docker.registryv2.flows.HttpFlowWithRetry.ContextWithRequest
import cromwell.docker.registryv2.flows.dockerhub.DockerHubFlow
import cromwell.docker.registryv2.flows.gcr.GoogleFlow
import cromwell.docker.registryv2.flows.quay.QuayFlow
import cromwell.engine.backend.{BackendSingletonCollection, CromwellBackends}
import cromwell.engine.io.{IoActor, IoActorProxy}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.AbortAllWorkflowsCommand
import cromwell.engine.workflow.lifecycle.execution.callcaching.{CallCache, CallCacheReadActor, CallCacheWriteActor}
import cromwell.engine.workflow.lifecycle.finalization.CopyWorkflowLogsActor
import cromwell.engine.workflow.tokens.{DynamicRateLimiter, JobExecutionTokenDispenserActor}
import cromwell.engine.workflow.workflowstore._
import cromwell.jobstore.{JobStore, JobStoreActor, SqlJobStore}
import cromwell.services.{EngineServicesStore, ServiceRegistryActor}
import cromwell.subworkflowstore.{SqlSubWorkflowStore, SubWorkflowStoreActor}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.concurrent.{Await, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

/**
  * An actor which serves as the lord protector for the rest of Cromwell, allowing us to have more fine grain
  * control on top level supervision, etc.
  *
  * For now there are only two (non-test) entry points into Cromwell - either using SingleWorkflowRunnerActor or CromwellServerActor.
  * This is intended to be extended by those entry points.
  *
  * If any of the actors created by CromwellRootActor fail to initialize the ActorSystem will die, which means that
  * Cromwell will fail to start in a bad state regardless of the entry point.
  * 
  * READ THIS: If you add a "system-level" actor here, make sure to consider what should be its
  * position in the shutdown process and modify CromwellShutdown accordingly.
  */
abstract class CromwellRootActor(gracefulShutdown: Boolean, abortJobsOnTerminate: Boolean)(implicit materializer: ActorMaterializer) extends Actor with ActorLogging with GracefulShutdownHelper {
  import CromwellRootActor._
  
  // Make sure the filesystems are initialized at startup
  val _ = CromwellFileSystems.instance

  private val logger = Logging(context.system, this)
  private val config = ConfigFactory.load()
  private implicit val system = context.system

  private val workflowHeartbeatConfig = WorkflowHeartbeatConfig(config)
  logger.info("Workflow heartbeat configuration:\n{}", workflowHeartbeatConfig)

  val serverMode: Boolean

  lazy val systemConfig = config.getConfig("system")
  lazy val serviceRegistryActor: ActorRef = context.actorOf(ServiceRegistryActor.props(config), "ServiceRegistryActor")
  lazy val numberOfWorkflowLogCopyWorkers = systemConfig.as[Option[Int]]("number-of-workflow-log-copy-workers").getOrElse(DefaultNumberOfWorkflowLogCopyWorkers)

  lazy val workflowStore: WorkflowStore = SqlWorkflowStore(EngineServicesStore.engineDatabaseInterface)
  lazy val workflowStoreCoordinatedWriteActor: ActorRef = context.actorOf(WorkflowStoreCoordinatedWriteActor.props(workflowStore))
  lazy val workflowStoreActor =
    context.actorOf(WorkflowStoreActor.props(
      workflowStoreDatabase = workflowStore,
      workflowStoreCoordinatedWriteActor = workflowStoreCoordinatedWriteActor,
      serviceRegistryActor = serviceRegistryActor,
      abortAllJobsOnTerminate = abortJobsOnTerminate,
      workflowHeartbeatConfig = workflowHeartbeatConfig),
      "WorkflowStoreActor")

  lazy val jobStore: JobStore = new SqlJobStore(EngineServicesStore.engineDatabaseInterface)
  lazy val jobStoreActor = context.actorOf(JobStoreActor.props(jobStore, serviceRegistryActor), "JobStoreActor")

  lazy val subWorkflowStore = new SqlSubWorkflowStore(EngineServicesStore.engineDatabaseInterface)
  lazy val subWorkflowStoreActor = context.actorOf(SubWorkflowStoreActor.props(subWorkflowStore), "SubWorkflowStoreActor")

  // Io Actor
  lazy val throttleElements = systemConfig.as[Option[Int]]("io.number-of-requests").getOrElse(100000)
  lazy val throttlePer = systemConfig.as[Option[FiniteDuration]]("io.per").getOrElse(100 seconds)
  lazy val nioParallelism = systemConfig.as[Option[Int]]("io.nio.parallelism").getOrElse(10)
  lazy val gcsParallelism = systemConfig.as[Option[Int]]("io.gcs.parallelism").getOrElse(10)
  lazy val ioThrottle = Throttle(throttleElements, throttlePer, throttleElements)
  lazy val ioActor = context.actorOf(IoActor.props(LoadConfig.IoQueueSize, nioParallelism, gcsParallelism, Option(ioThrottle), serviceRegistryActor), "IoActor")
  lazy val ioActorProxy = context.actorOf(IoActorProxy.props(ioActor), "IoProxy")

  lazy val workflowLogCopyRouter: ActorRef = context.actorOf(RoundRobinPool(numberOfWorkflowLogCopyWorkers)
    .withSupervisorStrategy(CopyWorkflowLogsActor.strategy)
    .props(CopyWorkflowLogsActor.props(serviceRegistryActor, ioActor)),
    "WorkflowLogCopyRouter")

  lazy val callCache: CallCache = new CallCache(EngineServicesStore.engineDatabaseInterface)

  lazy val numberOfCacheReadWorkers = config.getConfig("system").as[Option[Int]]("number-of-cache-read-workers").getOrElse(DefaultNumberOfCacheReadWorkers)
  lazy val callCacheReadActor = context.actorOf(RoundRobinPool(numberOfCacheReadWorkers)
    .props(CallCacheReadActor.props(callCache, serviceRegistryActor)),
    "CallCacheReadActor")

  lazy val callCacheWriteActor = context.actorOf(CallCacheWriteActor.props(callCache, serviceRegistryActor), "CallCacheWriteActor")

  // Docker Actor
  lazy val ioEc = context.system.dispatchers.lookup(Dispatcher.IoDispatcher)
  lazy val dockerConf = DockerConfiguration.instance
  // Sets the number of requests that the docker actor will accept before it starts backpressuring (modulo the number of in flight requests)
  lazy val dockerActorQueueSize = 500

  lazy val dockerHttpPool = Http().superPool[ContextWithRequest[DockerHashContext]]()
  lazy val googleFlow = new GoogleFlow(dockerHttpPool, dockerConf.gcrApiQueriesPer100Seconds)(ioEc, materializer, system.scheduler)
  lazy val dockerHubFlow = new DockerHubFlow(dockerHttpPool)(ioEc, materializer, system.scheduler)
  lazy val quayFlow = new QuayFlow(dockerHttpPool)(ioEc, materializer, system.scheduler)
  lazy val dockerCliFlow = new DockerCliFlow()(ioEc, system.scheduler)
  lazy val dockerFlows = dockerConf.method match {
    case DockerLocalLookup => Seq(dockerCliFlow)
    case DockerRemoteLookup => Seq(dockerHubFlow, googleFlow, quayFlow)
  }

  lazy val dockerHashActor = context.actorOf(DockerHashActor.props(dockerFlows, dockerActorQueueSize,
    dockerConf.cacheEntryTtl, dockerConf.cacheSize)(materializer), "DockerHashActor")

  lazy val backendSingletons = CromwellBackends.instance.get.backendLifecycleActorFactories map {
    case (name, factory) => name -> (factory.backendSingletonActorProps(serviceRegistryActor) map { context.actorOf(_, s"$name-Singleton") })
  }
  lazy val backendSingletonCollection = BackendSingletonCollection(backendSingletons)

  lazy val rate = DynamicRateLimiter.Rate(systemConfig.as[Int]("job-rate-control.jobs"), systemConfig.as[FiniteDuration]("job-rate-control.per"))

  lazy val jobExecutionTokenDispenserActor = context.actorOf(JobExecutionTokenDispenserActor.props(serviceRegistryActor, rate), "JobExecutionTokenDispenser")

  lazy val workflowManagerActor = context.actorOf(
    WorkflowManagerActor.props(
      workflowStore = workflowStoreActor,
      ioActor = ioActorProxy,
      serviceRegistryActor = serviceRegistryActor,
      workflowLogCopyRouter = workflowLogCopyRouter,
      jobStoreActor = jobStoreActor,
      subWorkflowStoreActor = subWorkflowStoreActor,
      callCacheReadActor = callCacheReadActor,
      callCacheWriteActor = callCacheWriteActor,
      dockerHashActor = dockerHashActor,
      jobTokenDispenserActor = jobExecutionTokenDispenserActor,
      backendSingletonCollection = backendSingletonCollection,
      serverMode = serverMode,
      workflowHeartbeatConfig = workflowHeartbeatConfig),
    "WorkflowManagerActor")

  if (gracefulShutdown) {
    // If abortJobsOnTerminate is true, aborting all workflows will be handled by the graceful shutdown process
    CromwellShutdown.registerShutdownTasks(
      cromwellId = workflowHeartbeatConfig.cromwellId,
      abortJobsOnTerminate,
      actorSystem = context.system,
      workflowManagerActor = workflowManagerActor,
      logCopyRouter = workflowLogCopyRouter,
      jobStoreActor = jobStoreActor,
      jobTokenDispenser = jobExecutionTokenDispenserActor,
      workflowStoreActor = workflowStoreActor,
      subWorkflowStoreActor = subWorkflowStoreActor,
      callCacheWriteActor = callCacheWriteActor,
      ioActor = ioActorProxy,
      dockerHashActor = dockerHashActor,
      serviceRegistryActor = serviceRegistryActor,
      materializer = materializer
    )
  } else if (abortJobsOnTerminate) {
    // If gracefulShutdown is false but abortJobsOnTerminate is true, set up a classic JVM shutdown hook
    sys.addShutdownHook {
      implicit val ec = context.system.dispatcher

      val abortFuture: Future[Unit] = for {
        // Give 30 seconds to the workflow store to switch all running workflows to aborting and shutdown. Should be more than enough
        _ <- gracefulStop(workflowStoreActor, 30.seconds, ShutdownCommand)
        // Once all workflows are "Aborting" in the workflow store, ask the WMA to effectively abort all of them
        _ <- gracefulStop(workflowManagerActor, AbortTimeout, AbortAllWorkflowsCommand)
      } yield ()

      Try(Await.result(abortFuture, AbortTimeout)) match {
        case Success(_) => logger.info("All workflows aborted")
        case Failure(f) => logger.error("Failed to abort workflows", f)
      }
    }
  }

  override def receive = {
    case message => logger.error(s"Unknown message received by CromwellRootActor: $message")
  }

  /**
    * Validate that all of the direct children actors were successfully created, otherwise error out the initialization
    * of Cromwell by passing a Throwable to the guardian.
    */
  override val supervisorStrategy = OneForOneStrategy() {
    case _: ActorInitializationException => Escalate
    case _: ActorRestartException => Restart
    case t => super.supervisorStrategy.decider.applyOrElse(t, (_: Any) => Escalate)
  }
}

object CromwellRootActor extends GracefulStopSupport {
  import net.ceedubs.ficus.Ficus._
  val AbortTimeout = ConfigFactory.load().as[FiniteDuration]("akka.coordinated-shutdown.phases.abort-all-workflows.timeout")
  val DefaultNumberOfWorkflowLogCopyWorkers = 10
  val DefaultCacheTTL = 20 minutes
  val DefaultNumberOfCacheReadWorkers = 25
}
