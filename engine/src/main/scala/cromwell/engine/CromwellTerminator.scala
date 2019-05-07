package cromwell.engine

import akka.Done
import akka.actor.CoordinatedShutdown

import scala.concurrent.Future

trait CromwellTerminator {
  def beginCromwellShutdown(reason: CoordinatedShutdown.Reason): Future[Done]
}
