package cromwell.engine.workflow.workflowstore

import akka.actor.{Actor, ActorLogging, ActorRef, Props, Timers}
import com.google.common.cache.CacheBuilder
import cromwell.core.{CacheConfig, WorkflowId}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.workflowstore.AbortRequestScanningActor.{AbortConfig, RunScan}
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.{FindWorkflowsWithAbortRequested, FindWorkflowsWithAbortRequestedFailure, FindWorkflowsWithAbortRequestedSuccess}

import scala.collection.JavaConverters._
import scala.concurrent.duration._


class AbortRequestScanningActor(workflowStoreActor: ActorRef,
                                workflowManagerActor: ActorRef,
                                abortConfig: AbortConfig,
                                workflowHeartbeatConfig: WorkflowHeartbeatConfig) extends Actor with Timers with ActorLogging {
  private val cache = {
    val cacheConfig = abortConfig.cacheConfig
    CacheBuilder.newBuilder()
      .concurrencyLevel(cacheConfig.concurrency)
      .expireAfterWrite(cacheConfig.ttl.length, cacheConfig.ttl.unit)
      .build[WorkflowId, java.lang.Boolean]()
  }

  self ! RunScan

  override def receive: Receive = {
    case RunScan =>
      workflowStoreActor ! FindWorkflowsWithAbortRequested(cromwellId = workflowHeartbeatConfig.cromwellId)
    case FindWorkflowsWithAbortRequestedSuccess(ids) =>
      val oldIds = cache.getAllPresent(ids.asJava).keySet().asScala
      // Diff this list to the cache.
      val newIds = ids.toSet -- oldIds
      if (newIds.nonEmpty) {
        // If there are new IDs then command the WorkflowManagerActor to abort them.
        workflowManagerActor ! WorkflowManagerActor.AbortWorkflowsCommand(newIds)
        newIds foreach { cache.put(_, true) }
      }
      setTimer()
    case FindWorkflowsWithAbortRequestedFailure(t) =>
      log.error(t, "Error searching for abort requests")
      // Don't dwell on it, restart that timer.
      setTimer()
  }

  private def setTimer(): Unit = timers.startSingleTimer(RunScan, RunScan, abortConfig.scanFrequency)
}

object AbortRequestScanningActor {
  case class AbortConfig(scanFrequency: FiniteDuration, cacheConfig: CacheConfig)

  case object RunScan

  def props(workflowStoreActor: ActorRef, workflowManagerActor: ActorRef, abortConfig: AbortConfig, workflowHeartbeatConfig: WorkflowHeartbeatConfig): Props =
    Props(new AbortRequestScanningActor(workflowStoreActor, workflowManagerActor, abortConfig, workflowHeartbeatConfig))
}
