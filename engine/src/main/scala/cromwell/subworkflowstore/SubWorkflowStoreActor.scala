package cromwell.subworkflowstore

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.ExecutionIndex._
import cromwell.core.{JobKey, MonitoringCompanionHelper, WorkflowId}
import cromwell.database.sql.tables.SubWorkflowStoreEntry
import cromwell.subworkflowstore.SubWorkflowStoreActor._

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

class SubWorkflowStoreActor(database: SubWorkflowStore) extends Actor with ActorLogging with MonitoringCompanionHelper {
  
  implicit val ec: ExecutionContext = context.dispatcher

  val subWorkflowStoreReceive: Receive = {
    case register: RegisterSubWorkflow => registerSubWorkflow(sender(), register)
    case query: QuerySubWorkflow => querySubWorkflow(sender(), query)
    case complete: WorkflowComplete => workflowComplete(sender(), complete)
  }
  
  override def receive = subWorkflowStoreReceive.orElse(monitoringReceive)
  
  private def registerSubWorkflow(replyTo: ActorRef, command: RegisterSubWorkflow) = {
    addWork()
    database.addSubWorkflowStoreEntry(
      command.rootWorkflowExecutionUuid.toString,
      command.parentWorkflowExecutionUuid.toString,
      command.jobKey.scope.fullyQualifiedName,
      command.jobKey.index.fromIndex,
      command.jobKey.attempt,
      command.subWorkflowExecutionUuid.toString
    ) onComplete { 
      case Success(_) =>
        replyTo ! SubWorkflowStoreRegisterSuccess(command) 
        removeWork()
      case Failure(ex) =>
        replyTo ! SubWorkflowStoreFailure(command, ex)
        removeWork()
    }
  }

  private def querySubWorkflow(replyTo: ActorRef, command: QuerySubWorkflow) = {
    val jobKey = command.jobKey
    database.querySubWorkflowStore(command.parentWorkflowExecutionUuid.toString, jobKey.scope.fullyQualifiedName, jobKey.index.fromIndex, jobKey.attempt) onComplete {
      case Success(Some(result)) => replyTo ! SubWorkflowFound(result)
      case Success(None) => replyTo ! SubWorkflowNotFound(command)
      case Failure(ex) => replyTo ! SubWorkflowStoreFailure(command, ex)
    }
  }

  private def workflowComplete(replyTo: ActorRef, command: WorkflowComplete) = {
    addWork()
    database.removeSubWorkflowStoreEntries(command.workflowExecutionUuid.toString) onComplete {
      case Success(_) => removeWork()
      case Failure(ex) => 
        replyTo ! SubWorkflowStoreFailure(command, ex)
        removeWork()
    }
  }
  
}

object SubWorkflowStoreActor {
  sealed trait SubWorkflowStoreActorCommand
  case class RegisterSubWorkflow(rootWorkflowExecutionUuid: WorkflowId, parentWorkflowExecutionUuid: WorkflowId, jobKey: JobKey, subWorkflowExecutionUuid: WorkflowId) extends SubWorkflowStoreActorCommand
  case class QuerySubWorkflow(parentWorkflowExecutionUuid: WorkflowId, jobKey: JobKey) extends SubWorkflowStoreActorCommand
  case class WorkflowComplete(workflowExecutionUuid: WorkflowId) extends SubWorkflowStoreActorCommand

  sealed trait SubWorkflowStoreActorResponse
  case class SubWorkflowStoreRegisterSuccess(command: RegisterSubWorkflow) extends SubWorkflowStoreActorResponse
  case class SubWorkflowFound(subWorkflowStoreEntry: SubWorkflowStoreEntry) extends SubWorkflowStoreActorResponse
  case class SubWorkflowNotFound(command: QuerySubWorkflow) extends SubWorkflowStoreActorResponse
  
  case class SubWorkflowStoreFailure(command: SubWorkflowStoreActorCommand, failure: Throwable) extends SubWorkflowStoreActorResponse
  
  def props(database: SubWorkflowStore) = Props(
    new SubWorkflowStoreActor(database)
  ).withDispatcher(EngineDispatcher)
}
