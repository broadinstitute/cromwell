package cromwell.docker

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.core.actor.RobustClientHelper
import cromwell.docker.DockerHashActor.DockerHashResponse

import scala.concurrent.duration.FiniteDuration

trait DockerClientHelper extends RobustClientHelper { this: Actor with ActorLogging =>
  
  protected def dockerHashingActor: ActorRef
  
  private [docker] def dockerResponseReceive: Receive = {
    case dockerResponse: DockerHashResponse if hasTimeout(dockerResponse.request) =>
      cancelTimeout(dockerResponse.request)
      receive.apply(dockerResponse)
    case (context: Any, dockerResponse: DockerHashResponse) if hasTimeout(context -> dockerResponse.request) =>
      cancelTimeout(context -> dockerResponse.request)
      receive.apply(context -> dockerResponse)
  }

  def sendDockerCommand(dockerHashRequest: DockerHashRequest, timeout: FiniteDuration = RobustClientHelper.DefaultRequestLostTimeout) = {
    robustSend(dockerHashRequest, dockerHashingActor, timeout)
  }

  def dockerReceive = robustReceive orElse dockerResponseReceive
}
