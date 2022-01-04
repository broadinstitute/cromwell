package cromwell.services.metadata.impl

import java.util.UUID
import java.util.concurrent.Callable

import com.google.common.cache.CacheBuilder
import cromwell.core.WorkflowId
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataString, MetadataValue}
import java.time.{Duration => JDuration}
import cromwell.services.metadata.impl.MetadataStatisticsRecorder._
import scala.concurrent.duration._
import scala.util.Try

object MetadataStatisticsRecorder {
  final case class HeavyMetadataAlert(workflowId: WorkflowId, count: Long)
  final case class WorkflowMetadataWriteStatistics(workflowId: WorkflowId, totalWrites: Long, lastLogged: Long, knownParent: Option[WorkflowId])
}

final class MetadataStatisticsRecorder(workflowCacheSize: Long = 100000L, // 100,000
                                       metadataAlertInterval: Long = 1L * 1000000, // 1 million
                                       bundleSubworkflowsIntoParents: Boolean = false
                                      ) {

  // Statistics for each workflow
  private val metadataWriteStatisticsCache = CacheBuilder.newBuilder()
    .concurrencyLevel(2)
    .expireAfterAccess(JDuration.ofSeconds(4.hours.toSeconds))
    .maximumSize(workflowCacheSize)
    .build[WorkflowId, WorkflowMetadataWriteStatistics]()

  def writeStatisticsLoader(workflowId: WorkflowId): Callable[WorkflowMetadataWriteStatistics] = () => WorkflowMetadataWriteStatistics(workflowId, 0L, 0L, None)

  def processEvents(putEvents: Iterable[MetadataEvent]): Vector[HeavyMetadataAlert] = {
    putEvents.groupBy(_.key.workflowId).toVector.flatMap { case (id, list) => processEventsForWorkflow(id, list)}
  }

  private def processEventsForWorkflow(workflowId: WorkflowId, events: Iterable[MetadataEvent]): Vector[HeavyMetadataAlert] = {
    val workflowWriteStats = metadataWriteStatisticsCache.get(workflowId, writeStatisticsLoader(workflowId))

    // Find a new parent record if one exists and update the statistics to record it:
    val parentallyUpdatedStatistics = {
      if (!bundleSubworkflowsIntoParents) workflowWriteStats
      else if (workflowWriteStats.knownParent.isDefined) workflowWriteStats
      else {
        val newParentId = events.collectFirst {
          case MetadataEvent(MetadataKey(_, None, "parentWorkflowId"), Some(MetadataValue(value, MetadataString)), _) => Try(WorkflowId(UUID.fromString(value))).toOption
        }.flatten
        workflowWriteStats.copy(knownParent = newParentId)
      }
    }

    updateStatisticsCacheAndGenerateAlerts(parentallyUpdatedStatistics, events.size.longValue)
  }

  private def updateStatisticsCacheAndGenerateAlerts(workflowWriteStats: WorkflowMetadataWriteStatistics, count: Long): Vector[HeavyMetadataAlert] = {
    val writesForWorkflow = workflowWriteStats.totalWrites + count

    val myAlerts = if (writesForWorkflow >= workflowWriteStats.lastLogged + metadataAlertInterval) {
      metadataWriteStatisticsCache.put(workflowWriteStats.workflowId, workflowWriteStats.copy(totalWrites = writesForWorkflow, lastLogged = writesForWorkflow))
      Vector(HeavyMetadataAlert(workflowWriteStats.workflowId, writesForWorkflow))
    }
    else {
      metadataWriteStatisticsCache.put(workflowWriteStats.workflowId, workflowWriteStats.copy(totalWrites = writesForWorkflow))
      Vector.empty
    }

    val parentalAlerts = workflowWriteStats.knownParent.toVector.flatMap { parentId =>
      val parentStatistics = metadataWriteStatisticsCache.get(parentId, writeStatisticsLoader(parentId))
      updateStatisticsCacheAndGenerateAlerts(parentStatistics, count)
    }

    myAlerts ++ parentalAlerts
  }

  // For testing/debugging only...:
  def statusString(): String = metadataWriteStatisticsCache.asMap().values().toString
}
