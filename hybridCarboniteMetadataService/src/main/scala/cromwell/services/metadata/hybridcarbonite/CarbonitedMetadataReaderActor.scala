package cromwell.services.metadata.hybridcarbonite

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.services.metadata.MetadataService.MetadataReadAction
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.FailedMetadataResponse

class CarbonitedMetadataReaderActor extends Actor with ActorLogging {

  override def receive: Receive = {
    case read: MetadataReadAction =>
      val error = new NotImplementedError("CarboniteWorkerActor")
      log.error(error, "Programmer Error! The CarbonitedMetadataReaderActor is not ready yet!")
      sender ! FailedMetadataResponse(read, error)
  }
}

object CarbonitedMetadataReaderActor {
  def props = Props(new CarbonitedMetadataReaderActor())
}
