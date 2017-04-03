package cromwell.server

import akka.actor.{ActorSystem, Terminated}
import akka.http.scaladsl.Http
import akka.stream.ActorMaterializer
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import org.slf4j.LoggerFactory

import scala.concurrent.Future

trait CromwellSystem {
  protected def systemName = "cromwell-system"
  protected def newActorSystem(): ActorSystem = ActorSystem(systemName)
  val conf = ConfigFactory.load()
  val logger = LoggerFactory.getLogger(getClass.getName)
  implicit final lazy val actorSystem = newActorSystem()
  implicit final lazy val materializer = ActorMaterializer()
  implicit private final lazy val ec = actorSystem.dispatcher

  def shutdownActorSystem(): Future[Terminated] = {
    Http().shutdownAllConnectionPools() flatMap { _ =>
      shutdownMaterializerAndActorSystem()
    } recoverWith {
      case _ => shutdownMaterializerAndActorSystem()
    }
  }
  
  private def shutdownMaterializerAndActorSystem() = {
    materializer.shutdown()
    actorSystem.terminate()
  }

  CromwellBackends.initBackends(BackendConfiguration.AllBackendEntries)
}
