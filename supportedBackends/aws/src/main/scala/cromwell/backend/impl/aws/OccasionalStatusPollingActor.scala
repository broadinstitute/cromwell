package cromwell.backend.impl.aws

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.backend.impl.aws.OccasionalStatusPollingActor._
import cromwell.backend.impl.aws.RunStatus.{Initializing, Running}
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.batch.model.ListJobsRequest

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

class OccasionalStatusPollingActor(configRegion: Option[Region]) extends Actor with ActorLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  val PollingInterval = 10.seconds

  lazy val client = {
    val builder = BatchClient.builder()
    configRegion.foreach(builder.region)
    builder.build
  }

  // Maps job ID to status
  var statuses: Map[String, RunStatus] = Map.empty

  override def receive = {
    case WhatsMyStatus(jobId) =>
      sender ! ThisWasYourStatus(statuses.get(jobId))

    case UpdateStatuses =>
      Future {
        updateStatuses()
      } onComplete {
        case Success(_) =>
          scheduleStatusUpdate(PollingInterval)
        case Failure(error) =>
          log.error(error, "Failed to update statuses. Will try again in 2 seconds")
          scheduleStatusUpdate(2.seconds)
      }
    case NotifyOfStatus(jobId, status) => statuses += jobId -> status
  }

  private def updateStatuses() = {

    def updateForStatus(awsStatusNames: Seq[String], runStatus: RunStatus) = {

      val valuesInStatus = awsStatusNames flatMap { awsStatusName =>
        val request = ListJobsRequest.builder()
          .jobStatus(awsStatusName)
          .jobQueue("GenomicsDefaultQueue-9d1df58e6cbf7dc")
          .build()


        client.listJobs(request).jobSummaryList().asScala
      }

      // Remove the old values and add the new values
      statuses = statuses.filterNot(_._2 == runStatus) ++ valuesInStatus.map(_.jobId() -> runStatus)
    }

    updateForStatus(List("SUBMITTED", "PENDING", "RUNNABLE"), Initializing)
    updateForStatus(List("STARTING", "RUNNING"), Running)
  }

  def scheduleStatusUpdate(in: FiniteDuration): Unit = {
    context.system.scheduler.scheduleOnce(in) { self ! UpdateStatuses }
    ()
  }

  scheduleStatusUpdate(PollingInterval)
}

object OccasionalStatusPollingActor {

  sealed trait OccasionalStatusPollingActorMessage

  case object UpdateStatuses extends OccasionalStatusPollingActorMessage
  case class StatusUpdates(newState: Map[String, RunStatus]) extends OccasionalStatusPollingActorMessage

  final case class NotifyOfStatus(jobId: String, runStatus: RunStatus) extends OccasionalStatusPollingActorMessage
  final case class WhatsMyStatus(jobId: String) extends OccasionalStatusPollingActorMessage
  final case class ThisWasYourStatus(status: Option[RunStatus]) extends OccasionalStatusPollingActorMessage

  def props(configRegion: Option[Region]) = Props(new OccasionalStatusPollingActor(configRegion))
}
