package cromwell.backend.impl.aws

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.backend.impl.aws.GoSlowJobSubmitActor.GoSlowSubmitActorMessage
import cromwell.backend.impl.aws.OccasionalStatusPollingActor.OccasionalStatusPollingActorMessage
import cromwell.core.Mailbox
import software.amazon.awssdk.regions.Region

class AwsBatchSingletonActor(configRegion: Option[Region]) extends Actor with ActorLogging {
  val awsOccasionalStatusPoller = context.actorOf(OccasionalStatusPollingActor.props(configRegion).withMailbox(Mailbox.PriorityMailbox), "AWSOccasionalStatusPoller")
  val awsGoSlowSubmitActor = context.actorOf(GoSlowJobSubmitActor.props(configRegion).withMailbox(Mailbox.PriorityMailbox), "AWSGoSlowSubmitter")

  override def receive = {

    case statusQuery: OccasionalStatusPollingActorMessage =>
      awsOccasionalStatusPoller.forward(statusQuery)
    case goSlowRequest: GoSlowSubmitActorMessage =>
      awsGoSlowSubmitActor.forward(goSlowRequest)
    case other =>
      log.error(s"Unknown message to AwsBatchSingletonActor: ${other.getClass.getSimpleName}. Dropping it.")
  }
}

object AwsBatchSingletonActor {
  def props(configRegion: Option[Region]) = Props(new AwsBatchSingletonActor(configRegion))
}
