package cromwell.webservice

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.Dispatcher
import cromwell.webservice.EngineStatsActor._
import scala.concurrent.duration._

/**
  * An imperfect collector of the number of workflows & jobs active in the system. Takes a list of WorkflowActor
  * ActorRefs and sends each one a message which will ultimately have the underlying WorkflowExecutionActor report
  * the number of active EJEAs.
  *
  * Because of the vagaries of timing, etc this is intended to give a rough idea of what's going on instead of
  * being a ground truth.
  */
final case class EngineStatsActor(workflowActors: List[ActorRef], replyTo: ActorRef, timeout: FiniteDuration) extends Actor with ActorLogging {
  implicit val ec = context.dispatcher

  private var jobCounts = Map.empty[ActorRef, Int]

  /*
    It's possible that WorkflowActors might disappear behind us and never manage to write us back.
    Instead of waiting longingly, watching a mailbox which might never receive some love instead wait
    a specified period of time and assume anything which was going to reply already has
  */
  val scheduledMsg = context.system.scheduler.scheduleOnce(timeout, self, ShutItDown)

  if (workflowActors.isEmpty) reportStats()
  else workflowActors foreach { _ ! JobCountQuery }

  override def receive = {
    case JobCount(count) =>
      jobCounts += (sender -> count)
      if (jobCounts.size == workflowActors.size) reportStats()
    case ShutItDown =>
      reportStats()
    case wompWomp => log.error("Unexpected message to EngineStatsActor: {}", wompWomp)
  }

  private def reportStats(): Unit = {
    replyTo ! EngineStats(jobCounts.size, jobCounts.values.sum)
    scheduledMsg.cancel()
    context stop self
  }
}

object EngineStatsActor {
  import scala.language.postfixOps

  def props(workflowActors: List[ActorRef], replyTo: ActorRef, timeout: FiniteDuration = MaxTimeToWait) = {
    Props(EngineStatsActor(workflowActors, replyTo, timeout)).withDispatcher(Dispatcher.ApiDispatcher)
  }

  sealed abstract class EngineStatsActorMessage
  private case object ShutItDown extends EngineStatsActorMessage
  case object JobCountQuery extends EngineStatsActorMessage
  final case class JobCount(count: Int) extends EngineStatsActorMessage
  val NoJobs = JobCount(0)

  final case class EngineStats(workflows: Int, jobs: Int)

  val MaxTimeToWait = 30 seconds
}
