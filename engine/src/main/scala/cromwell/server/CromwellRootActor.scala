package cromwell.server

import akka.actor.SupervisorStrategy.{Escalate, Restart}
import akka.actor.{Actor, ActorInitializationException, ActorRef, OneForOneStrategy}
import akka.event.Logging
import akka.http.scaladsl.Http
import akka.http.scaladsl.model.HttpRequest
import akka.routing.RoundRobinPool
import akka.stream.ActorMaterializer
import com.typesafe.config.ConfigFactory
import cromwell.core.Dispatcher
import cromwell.core.callcaching.docker.DockerHashActor
import cromwell.core.callcaching.docker.DockerHashActor.{DockerHashActorException, DockerHashContext}
import cromwell.core.callcaching.docker.registryv2.flows.dockerhub.DockerHubFlow
import cromwell.core.callcaching.docker.registryv2.flows.gcr.GoogleFlow
import cromwell.engine.backend.{BackendSingletonCollection, CromwellBackends}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.lifecycle.CopyWorkflowLogsActor
import cromwell.engine.workflow.lifecycle.execution.callcaching.{CallCache, CallCacheReadActor}
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor
import cromwell.engine.workflow.workflowstore.{SqlWorkflowStore, WorkflowStore, WorkflowStoreActor}
import cromwell.jobstore.{JobStore, JobStoreActor, SqlJobStore}
import cromwell.services.{ServiceRegistryActor, SingletonServicesStore}
import cromwell.subworkflowstore.{SqlSubWorkflowStore, SubWorkflowStoreActor}
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps

/**
  * An actor which serves as the lord protector for the rest of Cromwell, allowing us to have more fine grain
  * control on top level supervision, etc.
  *
  * For now there are only two (non-test) entry points into Cromwell - either using SingleWorkflowRunnerActor or CromwellServerActor.
  * This is intended to be extended by those entry points.
  *
  * If any of the actors created by CromwellRootActor fail to initialize the ActorSystem will die, which means that
  * Cromwell will fail to start in a bad state regardless of the entry point.
  */
 abstract class CromwellRootActor extends Actor {
  import CromwellRootActor._

  private val logger = Logging(context.system, this)
  private val config = ConfigFactory.load()
  private implicit val system = context.system
  private implicit val materializer = ActorMaterializer()
  
  val serverMode: Boolean

  lazy val serviceRegistryActor: ActorRef = context.actorOf(ServiceRegistryActor.props(config), "ServiceRegistryActor")
  lazy val numberOfWorkflowLogCopyWorkers = config.getConfig("system").as[Option[Int]]("number-of-workflow-log-copy-workers").getOrElse(DefaultNumberOfWorkflowLogCopyWorkers)

  lazy val workflowLogCopyRouter: ActorRef = context.actorOf(RoundRobinPool(numberOfWorkflowLogCopyWorkers)
      .withSupervisorStrategy(CopyWorkflowLogsActor.strategy)
      .props(CopyWorkflowLogsActor.props(serviceRegistryActor)),
      "WorkflowLogCopyRouter")

  lazy val workflowStore: WorkflowStore = SqlWorkflowStore(SingletonServicesStore.databaseInterface)
  lazy val workflowStoreActor = context.actorOf(WorkflowStoreActor.props(workflowStore, serviceRegistryActor), "WorkflowStoreActor")

  lazy val jobStore: JobStore = new SqlJobStore(SingletonServicesStore.databaseInterface)
  lazy val jobStoreActor = context.actorOf(JobStoreActor.props(jobStore), "JobStoreActor")

  lazy val subWorkflowStore = new SqlSubWorkflowStore(SingletonServicesStore.databaseInterface)
  lazy val subWorkflowStoreActor = context.actorOf(SubWorkflowStoreActor.props(subWorkflowStore), "SubWorkflowStoreActor")

  lazy val callCache: CallCache = new CallCache(SingletonServicesStore.databaseInterface)
  lazy val callCacheReadActor = context.actorOf(RoundRobinPool(25)
    .props(CallCacheReadActor.props(callCache)),
    "CallCacheReadActor")
  
  // Docker Actor
  lazy val ioEc = context.system.dispatchers.lookup(Dispatcher.IoDispatcher)
  lazy val gcrQueriesPer100Sec = config.getAs[Int]("docker.gcr-api-queries-per-100-seconds") getOrElse 1000
  lazy val dockerCacheEntryTTL = config.as[Option[FiniteDuration]]("docker.cache-entry-ttl").getOrElse(DefaultCacheTTL)
  lazy val dockerCacheSize = config.getAs[Long]("docker.cache-size") getOrElse 200L
  // Sets the number of requests that the docker actor will accept before it starts backpressuring (modulo the number of in flight requests)
  lazy val dockerActorQueueSize = 500
  
  lazy val dockerHttpPool = Http().superPool[(DockerHashContext, HttpRequest)]()
  lazy val googleFlow = new GoogleFlow(dockerHttpPool, gcrQueriesPer100Sec)(ioEc, materializer)
  lazy val dockerHubFlow = new DockerHubFlow(dockerHttpPool)(ioEc, materializer)
  lazy val dockerFlows = Seq(dockerHubFlow, googleFlow)
  lazy val dockerHashActor = context.actorOf(DockerHashActor.props(dockerFlows, dockerActorQueueSize, dockerCacheEntryTTL, dockerCacheSize)(materializer).withDispatcher(Dispatcher.IoDispatcher))

  lazy val backendSingletons = CromwellBackends.instance.get.backendLifecycleActorFactories map {
    case (name, factory) => name -> (factory.backendSingletonActorProps map context.actorOf)
  }
  lazy val backendSingletonCollection = BackendSingletonCollection(backendSingletons)

  lazy val jobExecutionTokenDispenserActor = context.actorOf(JobExecutionTokenDispenserActor.props)

  def abortJobsOnTerminate: Boolean

  lazy val workflowManagerActor = context.actorOf(
    WorkflowManagerActor.props(
      workflowStoreActor, serviceRegistryActor, workflowLogCopyRouter, jobStoreActor, subWorkflowStoreActor, callCacheReadActor,
      dockerHashActor, jobExecutionTokenDispenserActor, backendSingletonCollection, abortJobsOnTerminate, serverMode),
    "WorkflowManagerActor")

  override def receive = {
    case _ => logger.error("CromwellRootActor is receiving a message. It prefers to be left alone!")
  }

  /**
    * Validate that all of the direct children actors were successfully created, otherwise error out the initialization
    * of Cromwell by passing a Throwable to the guardian.
    */
  override val supervisorStrategy = OneForOneStrategy() {
    case actorInitializationException: ActorInitializationException => Escalate
    case dockerHash: DockerHashActorException => Restart
    case t => super.supervisorStrategy.decider.applyOrElse(t, (_: Any) => Escalate)
  }
}

object CromwellRootActor {
  val DefaultNumberOfWorkflowLogCopyWorkers = 10
  val DefaultCacheTTL = 20 minutes
}
