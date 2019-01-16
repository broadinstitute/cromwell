package cromwell.engine.workflow.tokens

import akka.actor.ActorRef
import cromwell.core.JobExecutionToken
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.TokenQueue.{DequeueResult, LeasedActor, TokenQueuePlaceholder}
import cromwell.engine.workflow.tokens.UnhoggableTokenPool.{ComeBackLater, Oink, TokenHoggingLease, TokensAvailable}
import io.github.andrebeat.pool.Lease
import spray.json.{JsArray, JsNumber, JsObject, JsString}

import scala.collection.immutable.Queue

object TokenQueue {
  case class DequeueResult(leasedActor: Option[LeasedActor], tokenQueue: TokenQueue)
  case class LeasedActor(queuePlaceholder: TokenQueuePlaceholder, lease: Lease[JobExecutionToken]) {
    def actor: ActorRef = queuePlaceholder.actor
  }
  def apply(tokenType: JobExecutionTokenType, logger: TokenEventLogger) = new TokenQueue(Map.empty, Vector.empty, logger, new UnhoggableTokenPool(tokenType))
  final case class TokenQueuePlaceholder(actor: ActorRef, hogGroup: String)
}

/**
  * A queue assigned to a pool.
  * Elements can be dequeued if the queue is not empty and there are tokens in the pool
  */
final case class TokenQueue(queues: Map[String, Queue[TokenQueuePlaceholder]],
                            queueOrder: Vector[String],
                            eventLogger: TokenEventLogger,
                            private [tokens] val pool: UnhoggableTokenPool) {
  val tokenType = pool.tokenType

  /**
    * Total number of entries in the queue
    */
  def size = queues.map(_._2.size).sum

  /**
    * Enqueues an actor
    *
    * @return the new token queue
    */
  def enqueue(placeholder: TokenQueuePlaceholder): TokenQueue = {
    if (queues.contains(placeholder.hogGroup)) {
      this.copy(
        queues = queues + (placeholder.hogGroup -> queues(placeholder.hogGroup).enqueue(placeholder))
      )
    } else {
      this.copy(
        queues = queues + (placeholder.hogGroup -> Queue[TokenQueuePlaceholder](placeholder)),
        queueOrder = queueOrder :+ placeholder.hogGroup
      )
    }
  }

  /**
    * Returns a dequeue'd actor if one exists and there's a token available for it
    * Returns an updated token queue based on this request (successful queuees get removed, hogs get sent to the back)
    */
  def dequeue: DequeueResult = recursingDequeue(queues, Vector.empty, queueOrder)

  private def recursingDequeue(queues: Map[String, Queue[TokenQueuePlaceholder]], queuesTried: Vector[String], queuesRemaining: Vector[String]): DequeueResult = {
    if (queuesRemaining.isEmpty) {
      DequeueResult(None, this)
    } else {
      val hogGroup = queuesRemaining.head
      val remainingHogGroups = queuesRemaining.tail
      val leaseTry = pool.tryAcquire(hogGroup)

      leaseTry match {
        case thl: TokenHoggingLease =>
          val oldQueue = queues(hogGroup)
          val (placeholder, newQueue) = oldQueue.dequeue
          val (newQueues, newQueueOrder) = if (newQueue.isEmpty) {
            (queues - hogGroup, remainingHogGroups ++ queuesTried)
          } else {
            (queues + (hogGroup -> newQueue), remainingHogGroups ++ queuesTried :+ hogGroup)
          }
          DequeueResult(Some(LeasedActor(placeholder, thl)), TokenQueue(newQueues, newQueueOrder, eventLogger, pool))
        case ComeBackLater =>
          // The pool is completely full right now, so there's no benefit trying the other hog groups:
          eventLogger.outOfTokens(tokenType.backend)
          DequeueResult(None, this)
        case Oink =>
          eventLogger.flagTokenHog(hogGroup)
          recursingDequeue(queues, queuesTried :+ hogGroup, remainingHogGroups)
      }
    }
  }

  def removeLostActor(lostActor: ActorRef): TokenQueue = this.copy(
    queues = queues.map { case (hogGroup, queue) =>
      hogGroup -> queue.filterNot(_.actor == lostActor)
    }
  )

  /**
    * Returns true if there's at least on element that can be dequeued, false otherwise.
    */
  def available: Boolean =
    queues.keys.exists { hg =>
      pool.available(hg) match {
        case TokensAvailable => true
        case ComeBackLater =>
          eventLogger.outOfTokens(tokenType.backend)
          false
        case Oink =>
          eventLogger.flagTokenHog(hg)
          false
      }
    }

  def statusPrint: JsObject = {

    val queueSizes = JsArray(queues.toVector.map { case (name, q) =>
      JsObject(Map(
        "name" -> JsString(name),
        "queue size" -> JsNumber(q.size)
      ))
    })

    JsObject(Map(
      "queue" -> queueSizes,
      "pool" -> pool.poolState
    ))
  }
}
