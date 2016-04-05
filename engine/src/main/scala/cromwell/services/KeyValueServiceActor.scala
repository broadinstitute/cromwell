package cromwell.services

import akka.actor.{Props, Actor}
import com.typesafe.config.Config
import cromwell.services.KeyValueServiceActor._

import scala.util.{Failure, Success, Try}

object KeyValueServiceActor {
  trait KeyValueServiceActorMessage
  case class Put(key: String, value: String) extends KeyValueServiceActorMessage
  case class PutSuccess(key: String, value: String) extends KeyValueServiceActorMessage
  case class PutFailure(key: String, value: String, failure: Throwable) extends KeyValueServiceActorMessage

  case class Get(key: String) extends KeyValueServiceActorMessage


  def props(serviceConfigPath: String, config: Config) = Props(KeyValueServiceActor(serviceConfigPath, config))
}

case class KeyValueServiceActor(serviceConfigPath: String, entireConfig: Config) extends Actor {


  def receive = {
    case Put(key, value) =>
      doPut(key, value) match {
        case Success(_) => sender ! PutSuccess (key, value)
        case Failure(ex: Throwable) => sender ! PutFailure(key, value, ex)
      }
  }

  private def doPut(key: String, value: String): Try[Unit] = ???
}
