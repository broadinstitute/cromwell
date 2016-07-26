package cromwell.server

import akka.actor.SupervisorStrategy.Escalate
import akka.actor.{Actor, ActorInitializationException, ActorRef, OneForOneStrategy, Props}
import akka.event.Logging
import akka.routing.RoundRobinPool
import com.typesafe.config.ConfigFactory
import cromwell.database.CromwellDatabase
import cromwell.engine.workflow.lifecycle.CopyWorkflowLogsActor
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.workflowstore.{SqlWorkflowStore, WorkflowStore, WorkflowStoreActor}
import cromwell.jobstore.JobStoreActor
import cromwell.services.ServiceRegistryActor
import lenthall.config.ScalaConfig.EnhancedScalaConfig

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

  lazy val serviceRegistryActor: ActorRef = context.actorOf(ServiceRegistryActor.props(config), "ServiceRegistryActor")
  lazy val numberOfWorkflowLogCopyWorkers = config.getConfig("system").getIntOr("number-of-workflow-log-copy-workers", default=DefaultNumberOfWorkflowLogCopyWorkers)

  lazy val workflowLogCopyRouter: ActorRef = context.actorOf(RoundRobinPool(numberOfWorkflowLogCopyWorkers)
      .withSupervisorStrategy(CopyWorkflowLogsActor.strategy)
      .props(CopyWorkflowLogsActor.props(serviceRegistryActor)),
      "WorkflowLogCopyRouter")

   lazy val workflowStore: WorkflowStore = SqlWorkflowStore(CromwellDatabase.databaseInterface)
   lazy val workflowStoreActor = context.actorOf(WorkflowStoreActor.props(workflowStore, serviceRegistryActor), "WorkflowStoreActor")

  lazy val jobStoreActor = context.actorOf(JobStoreActor.props, "JobStoreActor")

   lazy val workflowManagerActor = context.actorOf(WorkflowManagerActor.props(workflowStoreActor, serviceRegistryActor, workflowLogCopyRouter, jobStoreActor), "WorkflowManagerActor")

  override def receive = {
    case _ => logger.error("CromwellRootActor is receiving a message. It prefers to be left alone!")
  }

  /**
    * Validate that all of the direct children actors were successfully created, otherwise error out the initialization
    * of Cromwell by passing a Throwable to the guardian.
    */
  override val supervisorStrategy = OneForOneStrategy() {
    case aie: ActorInitializationException => throw new Throwable(s"Unable to create actor for ActorRef ${aie.getActor}", aie.getCause)
    case t => super.supervisorStrategy.decider.applyOrElse(t, (_: Any) => Escalate)
  }
}

object CromwellRootActor {
  val DefaultNumberOfWorkflowLogCopyWorkers = 10
}
