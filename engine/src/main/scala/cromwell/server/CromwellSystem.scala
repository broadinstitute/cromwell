package cromwell.server

import akka.Done
import akka.actor.{ActorSystem, CoordinatedShutdown, DeadLetter, Props, Terminated}
import akka.http.scaladsl.Http
import akka.stream.ActorMaterializer
import com.typesafe.config.Config
import cromwell.engine.CromwellTerminator
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import cromwell.languages.config.{CromwellLanguages, LanguageConfiguration}
import cromwell.services.{EngineServicesStore, MetadataServicesStore}

import scala.concurrent.Future

trait CromwellSystem extends CromwellTerminator {
  /*
  Initialize the service stores and their respective data access objects, if they haven't been already, before starting
  any services such as Metadata refresh and especially HTTP server binding. Some services like HTTP binding terminate
  the app when they failed to start, while the DAOs were still in the middle of non-atomic transactions such as
  database upgrades.

  This is the closest we have to:
  - https://github.com/broadinstitute/firecloud-orchestration/blob/3c9482b/DEVNOTES.md
  - https://github.com/broadinstitute/clio/blob/7253ec0/clio-server/src/main/scala/org/broadinstitute/clio/server/service/ServerService.scala#L58-L66
  - https://en.wikipedia.org/wiki/Data_access_object
   */
  EngineServicesStore.engineDatabaseInterface
  MetadataServicesStore.metadataDatabaseInterface

  protected def systemName = "cromwell-system"
  def config: Config

  protected def newActorSystem(): ActorSystem = {
    val system = ActorSystem(systemName, config)
    // This sets up our own dead letter listener that we can shut off during shutdown
    lazy val deadLetterListener = system.actorOf(Props(classOf[CromwellDeadLetterListener]), "DeadLetterListenerActor")
    system.eventStream.subscribe(deadLetterListener, classOf[DeadLetter])
    system
  }

  implicit final lazy val actorSystem = newActorSystem()
  implicit final lazy val materializer = ActorMaterializer()
  implicit private final lazy val ec = actorSystem.dispatcher

  override def beginCromwellShutdown(reason: CoordinatedShutdown.Reason): Future[Done] = {
    CromwellShutdown.instance(actorSystem).run(reason)
  }

  def shutdownActorSystem(): Future[Terminated] = {
    // If the actor system is already terminated it's already too late for a clean shutdown
    // Note: This does not protect again starting 2 shutdowns concurrently
    if (!actorSystem.whenTerminated.isCompleted) {
      Http().shutdownAllConnectionPools() flatMap { _ =>
        shutdownMaterializerAndActorSystem()
      } recoverWith {
        case _ =>
          // we still want to shutdown the materializer and actor system if shutdownAllConnectionPools failed
          shutdownMaterializerAndActorSystem()
      }
    } else actorSystem.whenTerminated
  }
  
  private def shutdownMaterializerAndActorSystem() = {
    materializer.shutdown()
    actorSystem.terminate()
  }

  CromwellBackends.initBackends(BackendConfiguration.AllBackendEntries)
  CromwellLanguages.initLanguages(LanguageConfiguration.AllLanguageEntries)
}
