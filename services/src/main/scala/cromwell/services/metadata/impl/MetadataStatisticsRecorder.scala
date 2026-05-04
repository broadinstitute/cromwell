package cromwell.services.metadata.impl

import java.util.UUID
import java.util.concurrent.Callable

import com.google.common.cache.CacheBuilder
import cromwell.core.WorkflowId
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataString, MetadataValue}
import java.time.{Duration => JDuration}
import net.ceedubs.ficus.Ficus._
import com.typesafe.config.Config
import cromwell.services.metadata.impl.MetadataStatisticsRecorder._

import scala.concurrent.duration._
import scala.util.Try

object MetadataStatisticsRecorder {
  sealed trait MetadataAlert {
    def workflowId: WorkflowId
    def count: Long
  }
  final case class HeavyMetadataAlert(workflowId: WorkflowId, count: Long) extends MetadataAlert
  final case class MaxMetadataAlert(workflowId: WorkflowId, count: Long) extends MetadataAlert
  final case class WorkflowMetadataWriteStatistics(workflowId: WorkflowId,
                                                   totalWrites: Long,
                                                   lastLogged: Long,
                                                   knownParent: Option[WorkflowId]
  )

  sealed trait MetadataStatisticsRecorderSettings
  case object MetadataStatisticsDisabled extends MetadataStatisticsRecorderSettings

  final case class MetadataStatisticsEnabled(workflowCacheSize: Long, metadataAlertInterval: Long, metadataLimit: Long)
      extends MetadataStatisticsRecorderSettings

  def apply(statisticsRecorderSettings: MetadataStatisticsRecorderSettings): MetadataStatisticsRecorder =
    statisticsRecorderSettings match {
      case MetadataStatisticsEnabled(cacheSize, interval, limit) =>
        new ActiveMetadataStatisticsRecorder(cacheSize, interval, limit)
      case MetadataStatisticsDisabled => new NoopMetadataStatisticsRecorder()
    }

  object MetadataStatisticsRecorderSettings {
    val defaultCacheSize = 20000L
    val defaultAlertInterval = 100000L
    val defaultLimit = 100000000L

    def apply(configSection: Option[Config]): MetadataStatisticsRecorderSettings =
      (configSection flatMap { conf: Config =>
        if (conf.as[Option[Boolean]]("enabled").forall(identity)) {
          val cacheSize: Long = conf.getOrElse("cache-size", defaultCacheSize)
          val metadataAlertInterval: Long = conf.getOrElse("metadata-row-alert-interval", defaultAlertInterval)
          val metadataLimit: Long = conf.getOrElse("metadata-row-limit", defaultAlertInterval)
          Option(MetadataStatisticsEnabled(cacheSize, metadataAlertInterval, metadataLimit))
        } else None

      }).getOrElse(MetadataStatisticsDisabled)
  }
}

sealed trait MetadataStatisticsRecorder {
  def processEventsAndGenerateAlerts(putEvents: Iterable[MetadataEvent]): Vector[MetadataAlert]
}

final class NoopMetadataStatisticsRecorder extends MetadataStatisticsRecorder {
  def processEventsAndGenerateAlerts(putEvents: Iterable[MetadataEvent]): Vector[MetadataAlert] = Vector.empty
}

final class ActiveMetadataStatisticsRecorder(workflowCacheSize: Long, metadataAlertInterval: Long, metadataLimit: Long)
    extends MetadataStatisticsRecorder {

  // Statistics for each workflow
  private val metadataWriteStatisticsCache = CacheBuilder
    .newBuilder()
    .expireAfterAccess(JDuration.ofSeconds(4.hours.toSeconds))
    .maximumSize(workflowCacheSize)
    .build[WorkflowId, WorkflowMetadataWriteStatistics]()

  def writeStatisticsLoader(workflowId: WorkflowId): Callable[WorkflowMetadataWriteStatistics] = () =>
    WorkflowMetadataWriteStatistics(workflowId, 0L, 0L, None)

  def processEventsAndGenerateAlerts(putEvents: Iterable[MetadataEvent]): Vector[MetadataAlert] =
    putEvents.groupBy(_.key.workflowId).toVector.flatMap { case (id, list) => processEventsForWorkflow(id, list) }

  private def processEventsForWorkflow(workflowId: WorkflowId,
                                       events: Iterable[MetadataEvent]
  ): Vector[MetadataAlert] = {
    val workflowWriteStats = metadataWriteStatisticsCache.get(workflowId, writeStatisticsLoader(workflowId))

    // Find a new parent record if one exists and update the statistics to record it:
    val parentallyUpdatedStatistics =
      if (workflowWriteStats.knownParent.isDefined) workflowWriteStats
      else {
        val newParentId = events.collectFirst {
          case MetadataEvent(MetadataKey(_, None, "parentWorkflowId"), Some(MetadataValue(value, MetadataString)), _) =>
            Try(WorkflowId(UUID.fromString(value))).toOption
        }.flatten
        workflowWriteStats.copy(knownParent = newParentId)
      }

    updateStatisticsCacheAndGenerateAlerts(parentallyUpdatedStatistics, events.size.longValue)
  }

  private def updateStatisticsCacheAndGenerateAlerts(workflowWriteStats: WorkflowMetadataWriteStatistics,
                                                     count: Long
  ): Vector[MetadataAlert] = {
    val writesForWorkflow = workflowWriteStats.totalWrites + count

    val myAlerts = if (writesForWorkflow >= workflowWriteStats.lastLogged + metadataAlertInterval) {
      metadataWriteStatisticsCache.put(
        workflowWriteStats.workflowId,
        workflowWriteStats.copy(totalWrites = writesForWorkflow, lastLogged = writesForWorkflow)
      )
      val heavyAlert = Vector(HeavyMetadataAlert(workflowWriteStats.workflowId, writesForWorkflow))

      // Check against limit once per interval.
      // Otherwise we would continuously spam the alert once its condition becomes true.
      // After we fail the workflow it should never reach its next interval.
      val maxAlert = if (writesForWorkflow > metadataLimit) {
        Vector(MaxMetadataAlert(workflowWriteStats.workflowId, writesForWorkflow))
      } else Vector.empty

      heavyAlert ++ maxAlert
    } else {
      metadataWriteStatisticsCache.put(workflowWriteStats.workflowId,
                                       workflowWriteStats.copy(totalWrites = writesForWorkflow)
      )
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
