package cromwell.server

import java.util.concurrent.atomic.AtomicBoolean

import akka.Done
import akka.actor._
import akka.http.scaladsl.Http
import akka.pattern.{AskTimeoutException, GracefulStopSupport}
import akka.routing.Broadcast
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cats.instances.future._
import cats.syntax.functor._
import cromwell.core.WorkflowId
import cromwell.core.WorkflowProcessingEvents.DescriptionEventValue.Released
import cromwell.engine.workflow.WorkflowManagerActor.{AbortAllWorkflowsCommand, PreventNewWorkflowsFromStarting}
import cromwell.engine.workflow.WorkflowProcessingEventPublishing
import cromwell.languages.util.ImportResolver.HttpResolver
import cromwell.services.{EngineServicesStore, MetadataServicesStore}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.slf4j.LoggerFactory

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

/**
  * Collection of methods and objects used to control Cromwell graceful shutdown process.
  */
object CromwellShutdown extends GracefulStopSupport {
  private val logger = LoggerFactory.getLogger("CromwellShutdown")

  // Includes DB writing actors, I/O Actor and DockerHashActor
  private val PhaseStopIoActivity = "stop-io-activity"
  // Shutdown phase allocated when "abort-jobs-on-terminate" is true to give time to the system to abort all workflows
  // This phase is at the same level as PhaseServiceRequestsDone in the dependency graph.
  private val PhaseAbortAllWorkflows = "abort-all-workflows"

  private val _shutdownInProgress: AtomicBoolean = new AtomicBoolean(false)
  private var _coordinatedShutdown: Option[CoordinatedShutdown] = None

  /**
    * Single instance of CoordinatedShutdown. Assumes only one ActorSystem.
    */
  def instance(implicit actorSystem: ActorSystem): CoordinatedShutdown = synchronized {
    _coordinatedShutdown match {
      case Some(v) => v
      case None =>
        val coordinatedShutdown = CoordinatedShutdown(actorSystem)
        _coordinatedShutdown = Option(coordinatedShutdown)
        coordinatedShutdown
    }
  }

  /**
    * Returns true of the coordinated shutdown process is in progress. False otherwise.
    */
  def shutdownInProgress(): Boolean = _shutdownInProgress.get()

  /**
    * Register a task to unbind from the port during the ServiceUnbind phase.
    */
  def registerUnbindTask(actorSystem: ActorSystem, serverBinding: Future[Http.ServerBinding]) = {
    instance(actorSystem).addTask(CoordinatedShutdown.PhaseServiceUnbind, "UnbindingServerPort") { () =>
      // At this point it's still safe to schedule work on the actor system's dispatcher
      implicit val ec = actorSystem.dispatcher
      for {
        binding <- serverBinding
        _ <- binding.unbind()
        _ = logger.info("Http server unbound.") 
      } yield Done
    }
  }

  /**
    * Register tasks on the coordinated shutdown instance allowing a controlled, ordered shutdown process
    * meant to prevent data loss and ensure consistency.
    * Calling this method will add a JVM shutdown hook.
    */
  def registerShutdownTasks(
                             cromwellId: String,
                             abortJobsOnTerminate: Boolean,
                             actorSystem: ActorSystem,
                             workflowManagerActor: ActorRef,
                             logCopyRouter: ActorRef,
                             jobTokenDispenser: ActorRef,
                             jobStoreActor: ActorRef,
                             workflowStoreActor: ActorRef,
                             subWorkflowStoreActor: ActorRef,
                             callCacheWriteActor: ActorRef,
                             ioActor: ActorRef,
                             dockerHashActor: ActorRef,
                             serviceRegistryActor: ActorRef,
                             materializer: ActorMaterializer
                           ): Unit = {

    val coordinatedShutdown = this.instance(actorSystem)

    def shutdownActor(actor: ActorRef,
                      phase: String,
                      message: AnyRef,
                      customTimeout: Option[FiniteDuration] = None)(implicit executionContext: ExecutionContext) = {
      coordinatedShutdown.addTask(phase, s"stop${actor.path.name.capitalize}") { () =>
        val timeout = coordinatedShutdown.timeout(phase)
        logger.info(s"Shutting down ${actor.path.name} - Timeout = ${timeout.toSeconds} seconds")

        val action = gracefulStop(actor, customTimeout.getOrElse(coordinatedShutdown.timeout(phase)), message)
        action onComplete {
          case Success(_) => logger.info(s"${actor.path.name} stopped")
          case Failure(_: AskTimeoutException) =>
            logger.error(s"Timed out trying to gracefully stop ${actor.path.name}. Forcefully stopping it.")
            actorSystem.stop(actor)
          case Failure(f) => logger.error(s"An error occurred trying to gracefully stop ${actor.path.name}.", f)
        }

        action map { _ => Done }
      }
    }

    implicit val ec = actorSystem.dispatcher

    // 1) Stop starting new workflows. This will stop the WMA from sending pull messages to the WorkflowStore
    coordinatedShutdown.addTask(CoordinatedShutdown.PhaseBeforeServiceUnbind, "stopWorkflowPulling") { () =>
      import akka.pattern.ask
      _shutdownInProgress.set(true)
      implicit val timeout = Timeout(coordinatedShutdown.timeout(CoordinatedShutdown.PhaseBeforeServiceUnbind))
      (workflowManagerActor ? PreventNewWorkflowsFromStarting) map { _ =>
        logger.info("Workflow polling stopped")
        Done
      }
    }

    /* 2) The socket is unbound from the port in the CoordinatedShutdown.PhaseServiceUnbind
      * See cromwell.engine.server.CromwellServer
     */

    /* 3) Finish processing all requests:
      *  - Release any WorkflowStore entries held by this Cromwell instance.
      *  - Publish workflow processing event metadata for the released workflows.
      *  - Stop the WorkflowStore: The port is not bound anymore so we can't have new submissions.
      *   Process what's left in the message queue and stop. 
      *   Note that it's possible that some submissions are still asynchronously being prepared at the 
      *   akka http API layer (CromwellApiService) to be sent to the WorkflowStore.
      *   Those submissions might be lost if the WorkflowStore shuts itself down when it's finished processing its current work.
      *   In that case the "ask" over in the CromwellAPIService will fail with a AskTimeoutException and should be handled appropriately.
      *   This process still ensures that no submission can make it to the database without a response being sent back to the client.
      *   
      *  - Stop WorkflowManagerActor: We've already stopped starting new workflows but the Running workflows are still
      *   going. This is tricky because all the actor hierarchy under the WMA can be in a variety of state combinations.
      *   Specifically there is an asynchronous gap in several cases between emission of messages towards engine level
      *   actors (job store, cache store, etc...) and emission of messages towards the metadata service for the same logical
      *   event (e.g: job complete).
      *   The current behavior upon restart however is to re-play the graph, skipping execution of completed jobs (determined by
      *   engine job store) but still re-submitting all related metadata events. This is likely sub-optimal, but is used here
      *   to simply stop the WMA (which will trigger all its descendants to be stopped recursively) without more coordination.
      *   Indeed even if the actor is stopped in between the above mentioned gap, metadata will be re-submitted anyway on restart,
      *   even for completed jobs.
      *   
      *  - Stop the LogCopyRouter: it can generate metadata events and must therefore be stopped before the service registry.
      *   Wrap the ShutdownCommand in a Broadcast message so the router forwards it to all its routees
      *   Use the ShutdownCommand because a PoisonPill could stop the routees in the middle of "transaction"
      *   with the IoActor. The routees handle the ShutdownCommand properly and shutdown only when they have
      *   no outstanding requests to the IoActor. When all routees are dead the router automatically stops itself.
      *   
      *  - Stop the job token dispenser: stop it before stopping WMA and its EJEA descendants because
      *  the dispenser is watching all EJEAs and would be flooded by Terminated messages otherwise  
    */
    coordinatedShutdown.addTask(CoordinatedShutdown.PhaseServiceRequestsDone, "releaseWorkflowStoreEntries") { () =>
      EngineServicesStore.engineDatabaseInterface.releaseWorkflowStoreEntries(cromwellId).map(count => {
        logger.info("{} workflows released by {}", count, cromwellId)
      }).as(Done)
    }
    coordinatedShutdown.addTask(CoordinatedShutdown.PhaseServiceRequestsDone, "publishMetadataForReleasedWorkflowStoreEntries") { () =>
      EngineServicesStore.engineDatabaseInterface.findWorkflows(cromwellId).map(ids =>
        ids foreach { id =>
          WorkflowProcessingEventPublishing.publish(WorkflowId.fromString(id), cromwellId, Released, serviceRegistryActor)
        }
      ).as(Done)
    }

    shutdownActor(workflowStoreActor, CoordinatedShutdown.PhaseServiceRequestsDone, ShutdownCommand)
    shutdownActor(logCopyRouter, CoordinatedShutdown.PhaseServiceRequestsDone, Broadcast(ShutdownCommand))
    shutdownActor(jobTokenDispenser, CoordinatedShutdown.PhaseServiceRequestsDone, ShutdownCommand)
    
    /*
      * Aborting is only a special case of shutdown. Instead of sending a PoisonPill, send a AbortAllWorkflowsCommand
      * Also attach this task to a special shutdown phase allowing for a longer timeout.
     */
    if (abortJobsOnTerminate) {
      val abortTimeout = coordinatedShutdown.timeout(PhaseAbortAllWorkflows)
      shutdownActor(workflowManagerActor, PhaseAbortAllWorkflows, AbortAllWorkflowsCommand, Option(abortTimeout))
    } else {
      // This is a pretty rough shutdown of the WMA, depending on which route we go (full let it crash, or more fine grained
      // state management) we might want to revisit this way of terminating the WMA and its descendants.
      shutdownActor(workflowManagerActor, CoordinatedShutdown.PhaseServiceRequestsDone, PoisonPill)
    }

    /* 4) Shutdown connection pools
      * This will close all akka http opened connection pools tied to the actor system.
      * The pools stop accepting new work but are given a chance to execute the work submitted prior to the shutdown call.
      * When this future returns, all outstanding connections to client will be terminated.
      * Note that this also includes connection pools like the one used to lookup docker hashes.
     */
    coordinatedShutdown.addTask(CoordinatedShutdown.PhaseServiceStop, "TerminatingConnections") { () =>
      Http(actorSystem).shutdownAllConnectionPools() as {
        logger.info("Connection pools shut down")
        Done
      }
    }

    /* 5) Stop system level actors that require writing to the database or I/O
      * - SubWorkflowStoreActor
      * - JobStoreActor
      * - CallCacheWriteActor
      * - ServiceRegistryActor
      * - DockerHashActor
      * - IoActor
    */
    List(subWorkflowStoreActor, jobStoreActor, callCacheWriteActor, serviceRegistryActor, dockerHashActor, ioActor) foreach {
      shutdownActor(_, PhaseStopIoActivity, ShutdownCommand)
    }

    // 6) Close database and stream materializer
    coordinatedShutdown.addTask(CoordinatedShutdown.PhaseBeforeActorSystemTerminate, "closeDatabase") { () =>
      EngineServicesStore.engineDatabaseInterface.close()
      MetadataServicesStore.metadataDatabaseInterface.close()
      logger.info("Database closed")
      Future.successful(Done)
    }
    coordinatedShutdown.addTask(CoordinatedShutdown.PhaseBeforeActorSystemTerminate, "shutdownMaterializer") { () =>
      materializer.shutdown()
      logger.info("Stream materializer shut down")
      Future.successful(Done)
    }

    // 7) Close out the backend used for WDL HTTP import resolution
    // http://sttp.readthedocs.io/en/latest/backends/start_stop.html
    coordinatedShutdown.addTask(CoordinatedShutdown.PhaseBeforeActorSystemTerminate, "wdlHttpImportResolverBackend") { () =>
      Future {
        HttpResolver.closeBackendIfNecessary()
        logger.info("WDL HTTP import resolver closed")
        Done
      }
    }
  }
}
