package cromwell.engine.workflow.tokens

import akka.actor.{Actor, ActorLogging, ActorRef, Props, Terminated}
import cromwell.core.JobExecutionToken
import JobExecutionToken._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor._
import cromwell.engine.workflow.tokens.TokenPool.TokenPoolPop

import scala.collection.immutable.Queue

class JobExecutionTokenDispenserActor extends Actor with ActorLogging {

  /**
    * Lazily created token pool. We only create a pool for a token type when we need it
    */
  var tokenPools: Map[JobExecutionTokenType, TokenPool] = Map.empty
  var tokenAssignments: Map[ActorRef, JobExecutionToken] = Map.empty

  override def receive: Actor.Receive = {
    case JobExecutionTokenRequest(tokenType) => sendTokenRequestResult(sender, tokenType)
    case JobExecutionTokenReturn(token) => unassign(sender, token)
    case Terminated(terminee) => onTerminate(terminee)
  }

  private def sendTokenRequestResult(sndr: ActorRef, tokenType: JobExecutionTokenType): Unit = {
    if (tokenAssignments.contains(sndr)) {
      sndr ! JobExecutionTokenDispensed(tokenAssignments(sndr))
    } else {
      context.watch(sndr)
      val updatedTokenPool = getTokenPool(tokenType).pop() match {
        case TokenPoolPop(newTokenPool, Some(token)) =>
          assignAndSendToken(sndr, token)
          newTokenPool
        case TokenPoolPop(sizedTokenPoolAndQueue: SizedTokenPoolAndActorQueue, None) =>
          val (poolWithActorEnqueued, positionInQueue) = sizedTokenPoolAndQueue.enqueue(sndr)
          sndr ! JobExecutionTokenDenied(positionInQueue)
          poolWithActorEnqueued
        case TokenPoolPop(someOtherTokenPool, None) =>
          //If this has happened, somebody's been playing around in this class and not covered this case:
          throw new RuntimeException(s"Unexpected token pool type didn't return a token: ${someOtherTokenPool.getClass.getSimpleName}")
      }

      tokenPools += tokenType -> updatedTokenPool
    }
  }

  private def getTokenPool(tokenType: JobExecutionTokenType): TokenPool = tokenPools.getOrElse(tokenType, createNewPool(tokenType))

  private def createNewPool(tokenType: JobExecutionTokenType): TokenPool = {
    val newPool = TokenPool(tokenType) match {
      case s: SizedTokenPool => SizedTokenPoolAndActorQueue(s, Queue.empty)
      case anythingElse => anythingElse
    }
    tokenPools += tokenType -> newPool
    newPool
  }

  private def assignAndSendToken(actor: ActorRef, token: JobExecutionToken) = {
    tokenAssignments += actor -> token
    actor ! JobExecutionTokenDispensed(token)
  }

  private def unassign(actor: ActorRef, token: JobExecutionToken): Unit = {
    if (tokenAssignments.contains(actor) && tokenAssignments(actor) == token) {
      tokenAssignments -= actor

      val pool = getTokenPool(token.jobExecutionTokenType) match {
        case SizedTokenPoolAndActorQueue(innerPool, queue) if queue.nonEmpty =>
          val (nextInLine, newQueue) = queue.dequeue
          assignAndSendToken(nextInLine, token)
          SizedTokenPoolAndActorQueue(innerPool, newQueue)
        case other =>
          other.push(token)
      }

      tokenPools += token.jobExecutionTokenType -> pool
      context.unwatch(actor)
      ()
    } else {
      log.error("Job execution token returned from incorrect actor: {}", token)
    }
  }

  private def onTerminate(terminee: ActorRef): Unit = {
    tokenAssignments.get(terminee) match {
      case Some(token) =>
        log.debug("Actor {} stopped without returning its Job Execution Token. Reclaiming it!", terminee)
        self.tell(msg = JobExecutionTokenReturn(token), sender = terminee)
      case None =>
        log.debug("Actor {} stopped while we were still watching it... but it doesn't have a token. Removing it from any queues if necessary", terminee)
        tokenPools = tokenPools map {
          case (tokenType, SizedTokenPoolAndActorQueue(pool, queue)) => tokenType -> SizedTokenPoolAndActorQueue(pool, queue.filterNot(_ == terminee))
          case (tokenType, other) => tokenType -> other
        }
    }
    context.unwatch(terminee)
    ()
  }
}

object JobExecutionTokenDispenserActor {

  def props = Props(new JobExecutionTokenDispenserActor).withDispatcher(EngineDispatcher)

  case class JobExecutionTokenRequest(jobExecutionTokenType: JobExecutionTokenType)
  case class JobExecutionTokenReturn(jobExecutionToken: JobExecutionToken)

  sealed trait JobExecutionTokenRequestResult
  case class JobExecutionTokenDispensed(jobExecutionToken: JobExecutionToken) extends JobExecutionTokenRequestResult
  case class JobExecutionTokenDenied(positionInQueue: Integer) extends JobExecutionTokenRequestResult

  case class SizedTokenPoolAndActorQueue(sizedPool: SizedTokenPool, queue: Queue[ActorRef]) extends TokenPool {
    override def currentLoans = sizedPool.currentLoans
    override def push(jobExecutionToken: JobExecutionToken) = SizedTokenPoolAndActorQueue(sizedPool.push(jobExecutionToken), queue)
    override def pop() = {
      val underlyingPop = sizedPool.pop()
      TokenPoolPop(SizedTokenPoolAndActorQueue(underlyingPop.newTokenPool.asInstanceOf[SizedTokenPool], queue), underlyingPop.poppedItem)
    }

    /**
      * Enqueues an actor (or just finds its current position)
      *
      * @return The actor's position in the queue
      */
    def enqueue(actor: ActorRef): (SizedTokenPoolAndActorQueue, Int) = {
      queue.indexOf(actor) match {
        case -1 =>
          val newQueue = queue :+ actor
          (SizedTokenPoolAndActorQueue(sizedPool, newQueue), newQueue.size - 1) // Convert from 1-indexed to 0-indexed
        case index => (this, index)
      }
    }
  }
}
