package cromwell.engine

import akka.Done
import akka.actor.CoordinatedShutdown

import scala.concurrent.Future

object MockCromwellTerminator extends CromwellTerminator {
  override def beginCromwellShutdown(notUsed: CoordinatedShutdown.Reason): Future[Done] = {
    Future.successful(Done)
  }
}
