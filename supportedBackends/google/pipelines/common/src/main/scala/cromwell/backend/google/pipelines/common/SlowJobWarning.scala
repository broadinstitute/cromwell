package cromwell.backend.google.pipelines.common

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging}
import cromwell.backend.google.pipelines.common.SlowJobWarning.{WarnAboutSlownessAfter, WarnAboutSlownessIfNecessary, WarningDetails}

import scala.concurrent.duration.FiniteDuration

trait SlowJobWarning { this: Actor with ActorLogging =>

  def slowJobWarningReceive: Actor.Receive = {
    case WarnAboutSlownessAfter(jobId, duration) =>
      alreadyWarned = false
      warningDetails = Some(WarningDetails(jobId, OffsetDateTime.now(), OffsetDateTime.now().plusSeconds(duration.toSeconds)))
    case WarnAboutSlownessIfNecessary => handleWarnMessage()
  }

  var warningDetails: Option[WarningDetails] = None
  var alreadyWarned: Boolean = false

  def warnAboutSlowJobIfNecessary(jobId: String) = {
    // Don't do anything here because we might need to update state.
    // Instead, send a message and handle this in the receive block.
    self ! WarnAboutSlownessIfNecessary
  }

  private def handleWarnMessage(): Unit = {
    if (!alreadyWarned) {
      warningDetails match {
        case Some(WarningDetails(jobId, startTime, warningTime)) if OffsetDateTime.now().isAfter(warningTime) =>
          log.warning(s"Job with ID '{}' has been running since {} and might not be making progress", jobId, startTime)
          alreadyWarned = true
        case _ => // Nothing to do
      }
    }
  }

}

object SlowJobWarning {
  final case class WarnAboutSlownessAfter(jobId: String, duration: FiniteDuration)
  case object WarnAboutSlownessIfNecessary

  final case class WarningDetails(jobId: String, startTime: OffsetDateTime, warningTime: OffsetDateTime)
}
