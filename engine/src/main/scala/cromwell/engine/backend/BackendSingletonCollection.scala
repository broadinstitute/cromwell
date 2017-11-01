package cromwell.engine.backend

import akka.actor.ActorRef

final case class BackendSingletonCollection(backendSingletonActors: Map[String, Option[ActorRef]])
