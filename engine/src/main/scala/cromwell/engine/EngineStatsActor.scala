package cromwell.engine

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.backend.{BackendJobDescriptorKey, JobExecutionMap}
import cromwell.core.ExecutionStatus.ExecutionStatus
import cromwell.core.{Dispatcher, JobKey, WorkflowId, WorkflowState}
import cromwell.engine.EngineStatsActor.{EngineStats, LogStats, StatsQuery}
import cromwell.engine.workflow.WorkflowActor.WorkflowStateTransition
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.WorkflowExecutionActorResponse
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionDiff
import cromwell.util.TerminalUtil

import scala.concurrent.duration._
import scala.language.postfixOps

final class EngineStatsActor(statusReportInterval: FiniteDuration) extends Actor with ActorLogging {
  implicit val ec = context.dispatcher

  private var workflows: Map[WorkflowId, WorkflowState] = Map.empty
  private var jobs: Map[JobKey, ExecutionStatus] = Map.empty
  
  private lazy val logHeader = List("Status", "Count")
  
  if (statusReportInterval != Duration.Zero) context.system.scheduler.scheduleOnce(statusReportInterval, self, LogStats)

  override def receive = {
    case StatsQuery => sender() ! EngineStats(summary(workflows), summary(jobs))
    case executionDiff: WorkflowExecutionDiff => updateJobs(executionDiff.backendJobChanges)
    case WorkflowStateTransition(id, state) => updateWorkflow(id, state)
    case weaResponse: WorkflowExecutionActorResponse => cleanupJobs(weaResponse.jobExecutionMap)
    case LogStats => 
      logStatus()
      context.system.scheduler.scheduleOnce(10 seconds, self, LogStats)
      ()
  }
  
  private def updateWorkflow(id: WorkflowId, workflowStatus: WorkflowState) = {
    if (workflowStatus.isTerminal) workflows -= id
    else workflows += (id -> workflowStatus)
  }
  
  private def updateJobs(statusChanges: Map[BackendJobDescriptorKey, ExecutionStatus]) = jobs ++= statusChanges
  
  private def summary[T](data: Map[_, T]): Map[String, Int] = {
    data.foldLeft(Map.empty[String, Int])((acc, v) => {
      val status = v._2
      val jobCount = acc.getOrElse(status.toString, 0) + 1
      acc updated (status.toString, jobCount)
    })
  }
  
  private def cleanupJobs(jobExecutionMap: JobExecutionMap) = {
    jobs = jobs filterNot {
      case (k, _) => jobExecutionMap.values.flatten.toList.contains(k)
    }
  }

  private def logStatus() = {
    val formattedJobs = {
      val jobsList = summary(jobs).toList map {
        case (status, count) => List(status, count.toString)
      }

      TerminalUtil.mdTable(jobsList, logHeader)
    }

    log.info(s"\n[Job Execution]\n$formattedJobs")
  }
}

object EngineStatsActor {
  import scala.language.postfixOps

  def props(statusReport: FiniteDuration = Duration.Zero) = {
    Props(new EngineStatsActor(statusReport)).withDispatcher(Dispatcher.ApiDispatcher)
  }

  sealed abstract class EngineStatsActorMessage
  case object StatsQuery extends EngineStatsActorMessage
  case object LogStats extends EngineStatsActorMessage

  final case class EngineStats(workflows: Map[String, Int], jobs: Map[String, Int])

  val MaxTimeToWait = 3 seconds
}
