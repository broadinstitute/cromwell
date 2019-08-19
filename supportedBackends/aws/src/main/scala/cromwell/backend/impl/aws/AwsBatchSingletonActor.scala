package cromwell.backend.impl.aws

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.backend.impl.aws.IntervalLimitedAwsJobSubmitActor.IntervalLimitedAwsJobSubmitActorMessage
import cromwell.backend.impl.aws.OccasionalStatusPollingActor.OccasionalStatusPollingActorMessage
import cromwell.cloudsupport.aws.auth.AwsAuthMode
import cromwell.core.Mailbox
import software.amazon.awssdk.regions.Region

class AwsBatchSingletonActor(configRegion: Option[Region], optAwsAuthMode: Option[AwsAuthMode] = None) extends Actor with ActorLogging {
  val awsOccasionalStatusPoller = context.actorOf(OccasionalStatusPollingActor.props(configRegion, optAwsAuthMode).withMailbox(Mailbox.PriorityMailbox), "OccasionalStatusPollingActor")
  val awsIntervalLimitedSubmitActor = context.actorOf(IntervalLimitedAwsJobSubmitActor.props(configRegion).withMailbox(Mailbox.PriorityMailbox), "IntervalLimitedAWSSubmitActor")

  override def receive = {

    case statusQuery: OccasionalStatusPollingActorMessage =>
      awsOccasionalStatusPoller.forward(statusQuery)
    case submitActorMessage: IntervalLimitedAwsJobSubmitActorMessage =>
      awsIntervalLimitedSubmitActor.forward(submitActorMessage)
    case other =>
      log.error("Unknown message to AwsBatchSingletonActor: {}. Dropping it.", other)
  }
}

object AwsBatchSingletonActor {
  def props(configRegion: Option[Region], optAwsAuthMode: Option[AwsAuthMode] = None) = Props(new AwsBatchSingletonActor(configRegion, optAwsAuthMode))
}
