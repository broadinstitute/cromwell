package cromwell.services

import akka.actor.{ActorRef, Actor}
import com.typesafe.config.Config
import ServiceRegistryActor._

object ServiceRegistryActor {

  /**
    * @return The list of services and paths to their config block.
     */
  def getServices: Map[String, String] = ??? // TODO: Read from config, e.g. ("KV"->"services.KV", "PubSub"->"services.PubSub")

  case class Send(serviceName: String, message: Any)
  case class SendFailure(serviceName: String)
}

class ServiceRegistryActor extends Actor {

  val globalConfig: Config = ???
  val services: Map[String, ActorRef] = getServices map { case (key,value) => (key, createServiceActor(value)) }

  def createServiceActor(configBlock: String): ActorRef = {
    // look in config block to get actor name
    // create actor with Props(configBlock, globalConfig)
    // Becomes supervisor of actor
    ???
  }

  def receive = {
    case Send(serviceName, message) =>
      services.get(serviceName) match {
        case Some(ref) => ref ! (message, sender)
        case None => sender ! SendFailure(serviceName)
      }
    case _ => ??? // Something sensible
  }
}
