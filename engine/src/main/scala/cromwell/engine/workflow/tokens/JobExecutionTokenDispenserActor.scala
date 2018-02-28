package cromwell.engine.workflow.tokens

import akka.actor.{Actor, ActorLogging, ActorRef, Props, Terminated, Timers}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.JobExecutionToken._
import cromwell.core.{ExecutionStatus, JobExecutionToken}
import cromwell.engine.instrumentation.JobInstrumentation
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor._
import cromwell.engine.workflow.tokens.TokenQueue.LeasedActor
import cromwell.services.instrumentation.CromwellInstrumentation._
import cromwell.services.instrumentation.CromwellInstrumentationScheduler
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import io.github.andrebeat.pool.Lease
import scala.concurrent.duration._

class JobExecutionTokenDispenserActor(override val serviceRegistryActor: ActorRef, rate: Rate) extends Actor with ActorLogging 
  with JobInstrumentation with CromwellInstrumentationScheduler with Timers {

  /**
    * Lazily created token queue. We only create a queue for a token type when we need it
    */
  var tokenQueues: Map[JobExecutionTokenType, TokenQueue] = Map.empty
  var tokenAssignments: Map[ActorRef, Lease[JobExecutionToken]] = Map.empty

  scheduleInstrumentation { 
    sendGaugeJob(ExecutionStatus.Running.toString, tokenAssignments.size.toLong)
    sendGaugeJob(ExecutionStatus.QueuedInCromwell.toString, tokenQueues.values.map(_.size).sum.toLong)
  }

  override def preStart() = {
    // Fixed for now, but soon to be modifiable 
    timers.startPeriodicTimer(TokensTimerKey, TokensAvailable(rate.n), rate.per)
    super.preStart()
  }

  override def receive: Actor.Receive = {
    case JobExecutionTokenRequest(tokenType) => enqueue(sender, tokenType)
    case JobExecutionTokenReturn => release(sender)
    case TokensAvailable(n) => distribute(n)
    case Terminated(terminee) => onTerminate(terminee)
    case ShutdownCommand => context stop self
  }

  private def enqueue(sndr: ActorRef, tokenType: JobExecutionTokenType): Unit = {
    if (tokenAssignments.contains(sndr)) {
      sndr ! JobExecutionTokenDispensed
    } else {
      context.watch(sndr)
      val updatedTokenQueue = getTokenQueue(tokenType).enqueue(sndr)
      tokenQueues += tokenType -> updatedTokenQueue
    }
  }

  private def getTokenQueue(tokenType: JobExecutionTokenType): TokenQueue = {
    tokenQueues.getOrElse(tokenType, createNewQueue(tokenType))
  }

  private def createNewQueue(tokenType: JobExecutionTokenType): TokenQueue = {
    val newQueue = TokenQueue(tokenType)
    tokenQueues += tokenType -> newQueue
    newQueue
  }

  private def distribute(n: Int) = if (tokenQueues.nonEmpty) {
    val iterator = RoundRobinQueueIterator(tokenQueues.values.toList)

    val nextTokens = iterator.take(n)
    
    nextTokens.foreach({
      case LeasedActor(actor, lease) if !tokenAssignments.contains(actor) =>
        tokenAssignments = tokenAssignments + (actor -> lease)
        incrementJob("Started")
        actor ! JobExecutionTokenDispensed
      // Only one token per actor, if you've already got one, you don't get this token !
      case LeasedActor(_, lease) => lease.release()
    })
    
    tokenQueues = iterator.updatedQueues.map(queue => queue.tokenType -> queue).toMap
  }

  private def release(actor: ActorRef): Unit = {
    tokenAssignments.get(actor) match {
      case Some(leaedToken) =>
        tokenAssignments -= actor
        leaedToken.release()
        context.unwatch(actor)
        ()
      case None =>  log.error("Job execution token returned from incorrect actor: {}", actor.path.name)
    }
  }

  private def onTerminate(terminee: ActorRef): Unit = {
    tokenAssignments.get(terminee) match {
      case Some(_) =>
        log.debug("Actor {} stopped without returning its Job Execution Token. Reclaiming it!", terminee)
        self.tell(msg = JobExecutionTokenReturn, sender = terminee)
      case None =>
        log.debug("Actor {} stopped while we were still watching it... but it doesn't have a token. Removing it from any queues if necessary", terminee)
        tokenQueues = tokenQueues map {
          case (tokenType, tokenQueue @ TokenQueue(queue, _)) => tokenType -> tokenQueue.copy(queue = queue.filterNot(_ == terminee))
        }
    }
    context.unwatch(terminee)
    ()
  }
}

object JobExecutionTokenDispenserActor {
  case class Rate(n: Int, per: FiniteDuration)
  case object TokensTimerKey
  case class TokensAvailable(n: Int)

  def props(serviceRegistryActor: ActorRef, rate: Rate) = Props(new JobExecutionTokenDispenserActor(serviceRegistryActor, rate)).withDispatcher(EngineDispatcher)

  case class JobExecutionTokenRequest(jobExecutionTokenType: JobExecutionTokenType)

  case object JobExecutionTokenReturn
  case object JobExecutionTokenDispensed
}
