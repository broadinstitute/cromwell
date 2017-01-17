package cromwell.io

import java.io.IOException

import akka.actor.SupervisorStrategy.{Escalate, Restart}
import akka.actor.{Actor, ActorInitializationException, OneForOneStrategy, Props}
import akka.routing.SmallestMailboxPool
import com.typesafe.config.Config
import cromwell.core.Dispatcher
import cromwell.services.io.IoActorCommand
import net.ceedubs.ficus.Ficus._

/**
  * Wrapper for the IO router actor.
  * This is a sub optimal level of indirection as messages could directly be directed to the router
  * instead of having to go through this actor.
  * The reason is that the Io Actor is under the umbrella of the service registry which instantiates its services generically using reflection,
  * which can't work for a router actor.
  * Was the router IoActor to become a first class citizen, this wrapper should be removed.
  */
class ServiceRouterActor(serviceConfig: Config, globalConfig: Config) extends Actor {
  private val routerSupervisorStrategy = OneForOneStrategy() {
    case e: ActorInitializationException => Escalate
    // If an IO exception escapes the routee for some reason - try to restart it
    case ioException: IOException => Restart
    case t => super.supervisorStrategy.decider.applyOrElse(t, (_: Any) => Escalate)
  }

  /*
   * Note: The router, its routees, as well as this actor all run on the IoDispatcher.
   */
  val ioDispatcher = Dispatcher.IoDispatcher

  // Might be worth experimenting different types of routers / resizers / dispatchers ?
  // See http://doc.akka.io/docs/akka/current/scala/routing.html
  // It could also be using the "config version" of akka deployment
  val parallelism = serviceConfig.as[Option[Int]]("IO.parallelism").getOrElse(20)
  val configRouter = SmallestMailboxPool(parallelism, supervisorStrategy = routerSupervisorStrategy)
  val router = context.actorOf(configRouter.props(NioServiceWorker.props().withDispatcher(ioDispatcher)), "iORouter")

  def receive = {
    case command: IoActorCommand[_] => router forward command
  }

  object ServiceRouterActor {
    def props(serviceConfig: Config, globalConfig: Config) = Props(new ServiceRouterActor(serviceConfig, globalConfig)).withDispatcher(ioDispatcher)
  }
}