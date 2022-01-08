package centaur

import centaur.CromwellTracker._
import com.typesafe.scalalogging.StrictLogging
import cromwell.api.model.WorkflowMetadata
import cromwell.core.WorkflowProcessingEvent
import cromwell.core.WorkflowProcessingEvents.DescriptionEventValue
import io.circe.{Decoder, Encoder}
import io.circe.generic.semiauto._
import io.circe.optics.JsonPath._
import io.circe.parser.parse
import io.circe.syntax._
import org.apache.commons.math3.stat.inference.ChiSquareTest

import scala.language.postfixOps


case class CromwellTracker(backendCount: Int, configuredSignificance: Double) extends StrictLogging {
  var counts: Map[String, Int] = Map()
  def track(metadata: WorkflowMetadata): Unit = {

    // Get all the workflow processing events for this workflow and pick out the "PickedUp" events. There can be multiple of
    // these for a workflow if Cromwell was restarted while it was running.
    val pickupEvents = parse(metadata.value) map {
      root.workflowProcessingEvents.each.as[WorkflowProcessingEvent].getAll(_)
    }

    pickupEvents match {
      case Left(exception) => throw exception
      case Right(events) =>
        events filter { _.description == DescriptionEventValue.PickedUp.value } foreach { e =>
          synchronized {
            counts = counts + (e.cromwellId -> (counts.getOrElse(e.cromwellId, 0) + 1))
            logger.info("added count for Cromwell id {}: {}", e.cromwellId, counts(e.cromwellId))
          }
        }
    }
  }

  def assertHoricromtality(): Unit = {
    counts foreach { case (id, count) => logger.info(s"horicromtal count {}: {}", id, count) }

    val runs = counts.values.sum
    val expected: Array[Double] = (1 to backendCount) map { _ => runs.toDouble / backendCount } toArray
    val actual: Array[Long] = counts.values map { _.toLong } toArray

    val observedSignificance = new ChiSquareTest().chiSquareTest(expected, actual)
    logger.info(f"configured/observed horicromtal significance levels: $configuredSignificance%.4f/$observedSignificance%.4f", configuredSignificance, observedSignificance)

    if (observedSignificance < configuredSignificance) {
      val message = f"Failed horicromtal check: observed significance level $observedSignificance%.4f, minimum of $configuredSignificance%.4f was required"
      throw new RuntimeException(message)
    }
  }

  override def toString: String = this.asInstanceOf[CromwellTracker].asJson.spaces2
}

object CromwellTracker {
  implicit lazy val cromwellTrackerEncoder: Encoder[CromwellTracker] = deriveEncoder
  implicit lazy val workflowProcessingEventEncoder: Encoder[WorkflowProcessingEvent] = deriveEncoder
  implicit lazy val workflowProcessingEventDecoder: Decoder[WorkflowProcessingEvent] = deriveDecoder
}
