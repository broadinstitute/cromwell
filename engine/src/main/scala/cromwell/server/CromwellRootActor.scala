package cromwell.server

import akka.actor.SupervisorStrategy.{Escalate, Restart}
import akka.actor.{Actor, ActorInitializationException, ActorLogging, ActorRef, OneForOneStrategy}
import akka.event.Logging
import akka.pattern.GracefulStopSupport
import akka.routing.RoundRobinPool
import akka.stream.ActorMaterializer
import com.typesafe.config.Config
import cromwell.core._
import cromwell.core.actor.StreamActorHelper.ActorRestartException
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.core.io.Throttle
import cromwell.core.io.Throttle._
import cromwell.docker.DockerInfoActor
import cromwell.docker.local.DockerCliFlow
import cromwell.engine.CromwellTerminator
import cromwell.engine.backend.{BackendSingletonCollection, CromwellBackends}
import cromwell.engine.io.{IoActor, IoActorProxy}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.AbortAllWorkflowsCommand
import cromwell.engine.workflow.lifecycle.execution.callcaching.{CallCache, CallCacheReadActor, CallCacheWriteActor}
import cromwell.engine.workflow.lifecycle.finalization.CopyWorkflowLogsActor
import cromwell.engine.workflow.tokens.{DynamicRateLimiter, JobExecutionTokenDispenserActor}
import cromwell.engine.workflow.workflowstore.AbortRequestScanningActor.AbortConfig
import cromwell.engine.workflow.workflowstore._
import cromwell.jobstore.{JobStore, JobStoreActor, SqlJobStore}
import cromwell.services.{EngineServicesStore, ServiceRegistryActor}
import cromwell.subworkflowstore.{SqlSubWorkflowStore, SubWorkflowStore, SubWorkflowStoreActor}
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
abstract class CromwellRootActor(terminator: CromwellTerminator,
                                 gracefulShutdown: Boolean,
                                 abortJobsOnTerminate: Boolean,
                                 val serverMode: Boolean,
                                 protected val config: Config)
                                (implicit materializer: ActorMaterializer)
  extends Actor with ActorLogging with GracefulShutdownHelper {

  import CromwellRootActor._
  
  // Make sure the filesystems are initialized at startup
  val _ = CromwellFileSystems.instance

  private val logger = Logging(context.system, this)

  private val workflowHeartbeatConfig = WorkflowHeartbeatConfig(config)
  logger.info("Workflow heartbeat configuration:\n{}", workflowHeartbeatConfig)

  lazy val systemConfig = config.getConfig("system")
  lazy val serviceRegistryActor: ActorRef = context.actorOf(ServiceRegistryActor.props(config), "ServiceRegistryActor")
  lazy val numberOfWorkflowLogCopyWorkers = systemConfig.as[Option[Int]]("number-of-workflow-log-copy-workers").getOrElse(DefaultNumberOfWorkflowLogCopyWorkers)

  lazy val workflowStore: WorkflowStore = SqlWorkflowStore(EngineServicesStore.engineDatabaseInterface)

  val workflowStoreAccess: WorkflowStoreAccess = {
    val coordinatedWorkflowStoreAccess = config.as[Option[Boolean]]("system.coordinated-workflow-store-access")
    coordinatedWorkflowStoreAccess match {
      case Some(false) => UncoordinatedWorkflowStoreAccess(workflowStore)
      case _ =>
        val coordinatedWorkflowStoreAccessActor: ActorRef = context.actorOf(
          WorkflowStoreCoordinatedAccessActor.props(workflowStore),
          "WorkflowStoreCoordinatedAccessActor"
        )
        CoordinatedWorkflowStoreAccess(coordinatedWorkflowStoreAccessActor)
    }
  }

  lazy val workflowStoreActor =
    context.actorOf(WorkflowStoreActor.props(
      workflowStoreDatabase = workflowStore,
      workflowStoreAccess = workflowStoreAccess,
      serviceRegistryActor = serviceRegistryActor,
      terminator = terminator,
      abortAllJobsOnTerminate = abortJobsOnTerminate,
      workflowHeartbeatConfig = workflowHeartbeatConfig),
      "WorkflowStoreActor")

  lazy val jobStore: JobStore = new SqlJobStore(EngineServicesStore.engineDatabaseInterface)
  lazy val jobStoreActor = context.actorOf(JobStoreActor.props(jobStore, serviceRegistryActor), "JobStoreActor")

  lazy val subWorkflowStore: SubWorkflowStore = new SqlSubWorkflowStore(EngineServicesStore.engineDatabaseInterface)
  lazy val subWorkflowStoreActor = context.actorOf(SubWorkflowStoreActor.props(subWorkflowStore), "SubWorkflowStoreActor")

  // Io Actor
  lazy val nioParallelism = systemConfig.as[Option[Int]]("io.nio.parallelism").getOrElse(10)
  lazy val gcsParallelism = systemConfig.as[Option[Int]]("io.gcs.parallelism").getOrElse(10)
  lazy val ioThrottle = systemConfig.getAs[Throttle]("io").getOrElse(Throttle(100000, 100.seconds, 100000))
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

  lazy val dockerCliFlow = new DockerCliFlow()(ioEc)
  lazy val dockerFlows = dockerConf.method match {
    case DockerLocalLookup => Seq(dockerCliFlow)
    case DockerRemoteLookup => DockerInfoActor.remoteRegistriesFromConfig(DockerConfiguration.dockerHashLookupConfig)
  }

  lazy val dockerHashActor = context.actorOf(DockerInfoActor.props(dockerFlows, dockerActorQueueSize,
    dockerConf.cacheEntryTtl, dockerConf.cacheSize), "DockerHashActor")

  lazy val backendSingletons = CromwellBackends.instance.get.backendLifecycleActorFactories map {
    case (name, factory) => name -> (factory.backendSingletonActorProps(serviceRegistryActor) map { context.actorOf(_, s"$name-Singleton") })
  }
  lazy val backendSingletonCollection = BackendSingletonCollection(backendSingletons)

  lazy val rate = DynamicRateLimiter.Rate(systemConfig.as[Int]("job-rate-control.jobs"), systemConfig.as[FiniteDuration]("job-rate-control.per"))
  lazy val tokenLogInterval: Option[FiniteDuration] = systemConfig.as[Option[Int]]("hog-safety.token-log-interval-seconds").map(_.seconds)

  lazy val jobExecutionTokenDispenserActor = context.actorOf(JobExecutionTokenDispenserActor.props(serviceRegistryActor, rate, tokenLogInterval), "JobExecutionTokenDispenser")

  lazy val workflowManagerActor = context.actorOf(
    WorkflowManagerActor.props(
      config = config,
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

  val abortRequestScanningActor = {
    val abortConfigBlock = config.as[Option[Config]]("system.abort")

    val abortCacheConfig = CacheConfig.config(caching = abortConfigBlock.flatMap { _.as[Option[Config]]("cache") },
                                              defaultConcurrency = 1,
                                              defaultSize = 100000L,
                                              defaultTtl = 20 minutes)

    val abortConfig = AbortConfig(
      scanFrequency = abortConfigBlock.flatMap { _.as[Option[FiniteDuration]]("scan-frequency") } getOrElse (30 seconds),
      cacheConfig = abortCacheConfig
    )

    context.actorOf(
      AbortRequestScanningActor.props(
        abortConfig = abortConfig,
        workflowStoreActor = workflowStoreActor,
        workflowManagerActor = workflowManagerActor,
        workflowHeartbeatConfig = workflowHeartbeatConfig
      ),
      "AbortRequestScanningActor"
    )
  }

  if (gracefulShutdown) {
    // If abortJobsOnTerminate is true, aborting all workflows will be handled by the graceful shutdown process
    CromwellShutdown.registerShutdownTasks(
      cromwellId = workflowHeartbeatConfig.cromwellId,
      abortJobsOnTerminate = abortJobsOnTerminate,
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
    val abortTimeout = config.as[FiniteDuration]("akka.coordinated-shutdown.phases.abort-all-workflows.timeout")
    sys.addShutdownHook {
      implicit val ec = context.system.dispatcher

      val abortFuture: Future[Unit] = for {
        // Give 30 seconds to the workflow store to switch all running workflows to aborting and shutdown. Should be more than enough
        _ <- gracefulStop(workflowStoreActor, 30.seconds, ShutdownCommand)
        // Once all workflows are "Aborting" in the workflow store, ask the WMA to effectively abort all of them
        _ <- gracefulStop(workflowManagerActor, abortTimeout, AbortAllWorkflowsCommand)
      } yield ()

      Try(Await.result(abortFuture, abortTimeout)) match {
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
  val DefaultNumberOfWorkflowLogCopyWorkers = 10
  val DefaultCacheTTL = 20 minutes
  val DefaultNumberOfCacheReadWorkers = 25
}
