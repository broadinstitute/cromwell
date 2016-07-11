package cromwell.engine.workflow

import java.time.OffsetDateTime

import akka.actor.{Actor, Props}
import akka.event.Logging
import cromwell.core.{WorkflowId, WorkflowSourceFiles}
import cromwell.database.obj.WorkflowMetadataKeys
import cromwell.engine.workflow.WorkflowStore.WorkflowToStart
import cromwell.engine.workflow.WorkflowStoreActor._
import cromwell.services.MetadataServiceActor.PutMetadataAction
import cromwell.services.{MetadataEvent, MetadataKey, MetadataValue, ServiceRegistryClient}

import scala.language.postfixOps
import scalaz.NonEmptyList

object WorkflowStore {
  sealed trait WorkflowState {
    def isStartable: Boolean = {
      this match {
        case x: StartableState => true
        case _ => false
      }
    }
  }

  case object Running extends WorkflowState

  sealed trait StartableState extends WorkflowState
  case object Submitted extends StartableState
  case object Restartable extends StartableState

  final case class SubmittedWorkflow(id: WorkflowId, sources: WorkflowSourceFiles, state: WorkflowState) {
    def toWorkflowToStart: WorkflowToStart = {
      state match {
        case r: StartableState => WorkflowToStart(id, sources, r)
        case _ => throw new IllegalArgumentException("This workflow is not currently in a startable state")
      }
    }
  }

  final case class WorkflowToStart(id: WorkflowId, sources: WorkflowSourceFiles, state: StartableState)
}

/**
  * Stores any submitted workflows which have not completed.
  *
  * This is not thread-safe, although it should not be being used in a multithreaded fashion. If that
  * changes, revisit.
  */
abstract class WorkflowStore {
  import WorkflowStore._

  var workflowStore = List.empty[SubmittedWorkflow]

  /**
    * Adds the requested WorkflowSourceFiles to the store and returns a WorkflowId for each one (in order)
    * for tracking purposes.
    */
  def add(sources: NonEmptyList[WorkflowSourceFiles]): NonEmptyList[WorkflowId] = {
    val submittedWorkflows = sources map { SubmittedWorkflow(WorkflowId.randomId(), _, Submitted) }
    workflowStore = workflowStore ++ submittedWorkflows.list
    submittedWorkflows map { _.id }
  }

  /**
    * Retrieves up to n workflows which have not already been pulled into the engine and sets their pickedUp
    * flag to true
    */
  def fetchRunnableWorkflows(n: Int): List[WorkflowToStart] = {
    val startableWorkflows = workflowStore filter { _.state.isStartable } take n
    val updatedWorkflows = startableWorkflows map { _.copy(state = Running) }
    workflowStore = (workflowStore diff startableWorkflows) ++ updatedWorkflows

    startableWorkflows map { _.toWorkflowToStart }
  }

  def remove(id: WorkflowId): Boolean = {
    val newWorkflowStore = workflowStore filterNot { _.id == id }

    if (newWorkflowStore == workflowStore) {
      false
    } else {
      workflowStore = newWorkflowStore
      true
    }
  }
}

class WorkflowStoreActor extends WorkflowStore with Actor with ServiceRegistryClient {
  /*
    WARNING: WorkflowStore is NOT thread safe. Unless that statement is no longer true do NOT use threads
    outside of the the single threaded Actor event loop
   */

  private val logger = Logging(context.system, this)

  override def receive = {
    case SubmitWorkflow(source) =>
      val id = add(NonEmptyList(source)).head
      registerIdWithMetadataService(id)
      sender ! WorkflowSubmittedToStore(id)
    case BatchSubmitWorkflows(sources) =>
      val ids = add(sources)
      ids foreach registerIdWithMetadataService
      sender ! WorkflowsBatchSubmittedToStore(ids)
    case FetchRunnableWorkflows(n) => sender ! newWorkflowMessage(n)
    case RemoveWorkflow(id) =>
      if (!remove(id)) logger.info(s"Attempted to remove ID $id from the WorkflowStore but it already exists!")
  }

  /**
    * Fetches at most n workflows, and builds the correct response message based on if there were any workflows or not
    */
  private def newWorkflowMessage(maxWorkflows: Int): WorkflowStoreActorResponse = {
    fetchRunnableWorkflows(maxWorkflows) match {
      case x :: xs => NewWorkflowsToStart(NonEmptyList.nel(x, xs))
      case _ => NoNewWorkflowsToStart
    }
  }

  /**
    * Takes the workflow id and sends it over to the metadata service w/ default empty values for inputs/outputs
    */
  private def registerIdWithMetadataService(id: WorkflowId): Unit = {
    val submissionEvents = List(
      MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionTime), MetadataValue(OffsetDateTime.now.toString)),
      MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Inputs)),
      MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Outputs))
    )

    serviceRegistryActor ! PutMetadataAction(submissionEvents)
  }
}

object WorkflowStoreActor {
  sealed trait WorkflowStoreActorCommand
  case class SubmitWorkflow(source: WorkflowSourceFiles) extends WorkflowStoreActorCommand
  case class BatchSubmitWorkflows(sources: NonEmptyList[WorkflowSourceFiles]) extends WorkflowStoreActorCommand
  case class FetchRunnableWorkflows(n: Int) extends WorkflowStoreActorCommand
  case class RemoveWorkflow(id: WorkflowId) extends WorkflowStoreActorCommand

  sealed trait WorkflowStoreActorResponse
  case class WorkflowSubmittedToStore(workflowId: WorkflowId) extends WorkflowStoreActorResponse
  case class WorkflowsBatchSubmittedToStore(workflowIds: NonEmptyList[WorkflowId]) extends WorkflowStoreActorResponse
  case object NoNewWorkflowsToStart extends WorkflowStoreActorResponse
  case class NewWorkflowsToStart(workflows: NonEmptyList[WorkflowToStart]) extends WorkflowStoreActorResponse

  def props() = Props(new WorkflowStoreActor)
}
