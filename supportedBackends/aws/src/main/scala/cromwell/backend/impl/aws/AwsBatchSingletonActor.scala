package cromwell.backend.impl.aws

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.backend.impl.aws.OccasionalStatusPollingActor.OccasionalStatusPollingActorMessage
import cromwell.core.Mailbox
import software.amazon.awssdk.regions.Region

class AwsBatchSingletonActor(configRegion: Option[Region]) extends Actor with ActorLogging {
  val awsOccasionalStatusPoller = context.actorOf(OccasionalStatusPollingActor.props(configRegion).withMailbox(Mailbox.PriorityMailbox), "AWSOccasionalStatusPoller")

  override def receive = {

    case statusQuery: OccasionalStatusPollingActorMessage =>
      awsOccasionalStatusPoller.forward(statusQuery)
    case other =>
      log.error(s"Unknown message to AwsBatchSingletonActor: ${other.getClass.getSimpleName}. Dropping it.")
  }
}

object AwsBatchSingletonActor {
  def props(configRegion: Option[Region]) = Props(new AwsBatchSingletonActor(configRegion))
}
