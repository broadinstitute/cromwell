package cromwell.engine.db

import akka.actor.{ActorRef, FSM, PoisonPill, Props}
import cromwell.engine.CromwellFSM
import cromwell.engine.CromwellFSM.ActorFailure
import cromwell.engine.db.DBQueryActor.DBQueryState
import cromwell.logging.WorkflowLogger

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

object DBQueryActor {
  sealed trait DBQueryState
  case object Idle extends DBQueryState
  case object Executing extends DBQueryState
  case object Done extends DBQueryState

  sealed trait DBQueryMessage
  case class UnitQueries(queries: List[Future[Unit]])
  case class YieldingQuery[T](query: Future[T])
  case class SuccessfulQuery[T](message: Any, result: T)
  case class FailedQuery(failureContext: String, failure: Throwable, message: Any) extends ActorFailure

  def props(replyTo: ActorRef, logger: WorkflowLogger, senderMessage: Any)(implicit ec: ExecutionContext) = {
    Props(new DBQueryActor(replyTo, logger, senderMessage))
  }
}

class DBQueryActor(replyTo: ActorRef, val logger: WorkflowLogger, senderMessage: Any)(implicit ec: ExecutionContext) extends CromwellFSM[DBQueryState, Unit] {
  import cromwell.engine.db.DBQueryActor._

  startWith(Idle, ())

  when(Idle) {
    case Event(UnitQueries(queries), _) =>
      Future.sequence(queries) onComplete {
        case Success(_) =>
          replyTo ! SuccessfulQuery(senderMessage, ())
          self ! PoisonPill
        case Failure(f) => fail(FailedQuery("Failed to perform DB operation", f, senderMessage))
      }

      goto(Executing)

    case Event(YieldingQuery(query), _) =>
      query onComplete {
        case Success(result) =>
          replyTo ! SuccessfulQuery(senderMessage, result)
          self ! PoisonPill
        case Failure(f) => fail(FailedQuery("Failed to perform DB operation", f, senderMessage))
      }

      goto(Executing)
  }

  when(Executing) { FSM.NullFunction }
  when(Done) { FSM.NullFunction }

  override def failureState: DBQueryState = Done
  override def sendFailureTo: Option[ActorRef] = Option(replyTo)
}
