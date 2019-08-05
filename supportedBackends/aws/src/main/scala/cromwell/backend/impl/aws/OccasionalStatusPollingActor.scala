package cromwell.backend.impl.aws

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.backend.impl.aws.OccasionalStatusPollingActor._
import cromwell.backend.impl.aws.RunStatus.{Initializing, Running}
import cromwell.cloudsupport.aws.auth.AwsAuthMode
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.batch.model.ListJobsRequest

import scala.annotation.tailrec
import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

/**
  * Works by polling AWS for the statuses of jobs in our queues of interest on a fixed schedule.
  * When queried for status, returns the currently known status, rather than going to AWS again.
  *
  *
  * Known Limitations (for they are legion):
  *
  *   - No IO retry in case a request to AWS fails (we just give up and don't update *anything* for this cycle)
  *   - The requests to AWS are also blocking IO
  *   - The 'queuesToMonitor' get added to but never cleared out (so on a multi-tenant system, will grow indefinitely)
  *   - We don't track completed jobs - so when a job completes the caller will get a None, and have to fall back to an AWS query anyway.
  */
class OccasionalStatusPollingActor(configRegion: Option[Region], optAwsAuthMode: Option[AwsAuthMode] = None) extends Actor with ActorLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  val PollingInterval = 10.seconds

  lazy val client = {
    val builder = BatchClient.builder()
    optAwsAuthMode.foreach { awsAuthMode =>
      builder.credentialsProvider(StaticCredentialsProvider.create(awsAuthMode.credential(_ => "")))
    }
    configRegion.foreach(builder.region)
    builder.build
  }

  // Maps job ID to status
  var statuses: Map[String, RunStatus] = Map.empty

  var queuesToMonitor: Set[String] = Set.empty

  override def receive = {
    case WhatsMyStatus(queueArn, jobId) =>
      queuesToMonitor += queueArn // Set addition so expectation is a no-op almost every time
      sender ! NotifyOfStatus(queueArn, jobId, statuses.get(jobId))

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
    case NotifyOfStatus(queueArn, jobId, Some(status)) =>
      statuses += jobId -> status
      queuesToMonitor += queueArn // Set addition so expectation is a no-op almost every time
    case NotifyOfStatus(_, jobId, None) =>
      log.error("Programmer Error: OccasionalStatusPollerActor was given an empty status update for {}. It was probably intended to be filled.", jobId)
  }

  private def updateStatuses() = {

    final case class PageAccumulation(nextPageToken: String, currentList: Vector[String])

    @tailrec
    def findJobsInStatus(awsStatusName: String, queueName: String, pageAccumulation: Option[PageAccumulation]): Vector[String] = {

      val requestBuilder = ListJobsRequest.builder()
        .jobStatus(awsStatusName)
        .jobQueue(queueName)
          .maxResults(100)

      pageAccumulation.foreach { pa => requestBuilder.nextToken(pa.nextPageToken) }

      val request = requestBuilder.build()

      val response = client.listJobs(request)
      val jobIds = response.jobSummaryList().asScala.map(_.jobId()).toVector

      val allKnownJobIds = jobIds ++ pageAccumulation.toVector.flatMap(_.currentList)

      Option(response.nextToken()) match {
        case Some(nextPageToken) =>
          findJobsInStatus(awsStatusName, queueName, Option(PageAccumulation(nextPageToken, allKnownJobIds)))
        case None => allKnownJobIds
      }
    }

    def updateForStatusNames(awsStatusNamesToCheck: Seq[String], mapToRunStatus: RunStatus): Unit = {
      Try {
        val jobIdsInStatus = for {
          queueName <- queuesToMonitor
          awsStatusName <- awsStatusNamesToCheck
          jobId <- findJobsInStatus(awsStatusName, queueName, pageAccumulation = None)
        } yield jobId

        // Remove the old values and add the new values

        statuses = statuses.filterNot(_._2 == mapToRunStatus) ++ jobIdsInStatus.map(_ -> mapToRunStatus)
      } recover {
        case e =>
          log.error(e, s"Failure fetching statuses for AWS jobs in $mapToRunStatus. No updates will occur.")
      }
      ()
    }

      updateForStatusNames(List("SUBMITTED", "PENDING", "RUNNABLE"), Initializing)
      updateForStatusNames(List("STARTING", "RUNNING"), Running)
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

  final case class NotifyOfStatus(queueArn: String, jobId: String, runStatus: Option[RunStatus]) extends OccasionalStatusPollingActorMessage
  final case class WhatsMyStatus(queueArn: String, jobId: String) extends OccasionalStatusPollingActorMessage

  def props(configRegion: Option[Region], optAwsAuthMode: Option[AwsAuthMode] = None) = Props(new OccasionalStatusPollingActor(configRegion, optAwsAuthMode))
}
