
package cromwell.server

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import org.slf4j.LoggerFactory

trait CromwellSystem {
  protected def systemName = "cromwell-system"
  protected def newActorSystem(): ActorSystem = ActorSystem(systemName)
  val conf = ConfigFactory.load()
  val logger = LoggerFactory.getLogger(getClass.getName)
  implicit final lazy val actorSystem = newActorSystem()

  def shutdownActorSystem(): Unit = actorSystem.shutdown()

  CromwellBackends.initBackends(BackendConfiguration.AllBackendEntries, actorSystem)
}