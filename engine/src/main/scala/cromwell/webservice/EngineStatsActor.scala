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
   * FIXME 
   * Because of sub workflows there is currently no reliable way to know if we received responses from all running WEAs. 
   * For now, we always wait for the timeout duration before responding to give a chance to all WEAs to respond (even nested ones).
   * This could be improved by having WEAs wait for their sub WEAs before sending back the response.
  */
  val scheduledMsg = context.system.scheduler.scheduleOnce(timeout, self, ShutItDown)

  if (workflowActors.isEmpty) reportStats()
  else workflowActors foreach { _ ! JobCountQuery }

  override def receive = {
    case JobCount(count) =>
      jobCounts += (sender -> count)
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

  val MaxTimeToWait = 3 seconds
}
