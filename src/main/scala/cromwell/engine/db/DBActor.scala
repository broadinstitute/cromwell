package cromwell.engine.db

import akka.actor.{Actor, ActorRef, Props}
import cromwell.engine.ExecutionStatus._
import cromwell.engine._
import cromwell.engine.backend.BackendCall
import cromwell.engine.db.DBActor.{DBMessage, PersistCallCompleteData, PersistCallStatusTerminal}
import cromwell.engine.db.DBQueryActor.UnitQueries
import cromwell.engine.db.DataAccess._
import cromwell.engine.workflow.OutputKey
import cromwell.logging.WorkflowLogger
import cromwell.server.GlobalActorSystem

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.language.implicitConversions

object DBActor {
  sealed trait DBMessage {
    def logger: WorkflowLogger
  }
  case class PersistCallCompleteData(workflowId: WorkflowId,
                                     key: OutputKey,
                                     outputs: CallOutputs,
                                     events: Seq[ExecutionEventEntry],
                                     logger: WorkflowLogger) extends DBMessage
  case class PersistCallStatusTerminal(workflowId: WorkflowId,
                                       key: ExecutionDatabaseKey,
                                       status: ExecutionStatus,
                                       scriptReturnCode: Option[Int],
                                       hash: Option[ExecutionHash],
                                       resultsClonedFrom: Option[BackendCall],
                                       logger: WorkflowLogger) extends DBMessage

  val instance = GlobalActorSystem.globalActorSystem.actorOf(Props(new DBActor))
}

class DBActor extends Actor with CromwellActor {

  implicit def queryToList(query: Future[Unit]): List[Future[Unit]] = List(query)

  private def queryActor(replyTo: ActorRef, message: DBMessage) = {
    context.actorOf(DBQueryActor.props(replyTo, message.logger, message))
  }

  override def receive: Receive = {
    case m: PersistCallCompleteData =>
      val outputsFuture = globalDataAccess.setOutputs(m.workflowId, m.key, m.outputs, m.key.scope.rootWorkflow.outputs)
      val eventsFuture = globalDataAccess.setExecutionEvents(m.workflowId, m.key.scope.fullyQualifiedName, m.key.index, m.key.attempt, m.events)
      queryActor(sender(), m) ! UnitQueries(List(outputsFuture, eventsFuture))
    case m: PersistCallStatusTerminal =>
      val statusFuture = globalDataAccess.setTerminalStatus(m.workflowId, m.key, m.status, m.scriptReturnCode, m.hash, m.resultsClonedFrom)
      queryActor(sender(), m) ! UnitQueries(statusFuture)
  }

}
