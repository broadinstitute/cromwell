package cromwell.server

import akka.actor.{ActorSystem, Terminated}
import akka.http.scaladsl.Http
import akka.stream.ActorMaterializer
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import cromwell.services.SingletonServicesStore

import scala.concurrent.Future

trait CromwellSystem {
  /*
  Initialize the services-store (aka a collection of DAOs), if they haven't been already, before starting any services
  such as Metadata refresh and especially HTTP server binding. Some services like HTTP binding terminate the app when
  they failed to start, while the DAOs were still in the middle of non-atomic transaction such as a database upgrade.

  This is the closest we have to:
  - https://github.com/broadinstitute/firecloud-orchestration/blob/3c9482b/DEVNOTES.md
  - https://github.com/broadinstitute/clio/blob/7253ec0/clio-server/src/main/scala/org/broadinstitute/clio/server/service/ServerService.scala#L58-L66
  - https://en.wikipedia.org/wiki/Data_access_object
   */
  SingletonServicesStore

  protected def systemName = "cromwell-system"
  val conf = ConfigFactory.load()
  
  protected def newActorSystem(): ActorSystem = ActorSystem(systemName, conf)
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
