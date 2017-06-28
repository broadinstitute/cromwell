package cromwell.server

import akka.actor.{ActorContext, Props}
import akka.http.scaladsl.Http
import akka.http.scaladsl.server.Route
import akka.http.scaladsl.server.Directives._

import akka.stream.ActorMaterializer
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.webservice.{CromwellApiService, SwaggerService}
import cromwell.webservice.WrappedRoute._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.concurrent.{Await, Future}
import scala.util.{Failure, Success}

// Note that as per the language specification, this is instantiated lazily and only used when necessary (i.e. server mode)
object CromwellServer {
  def run(cromwellSystem: CromwellSystem): Future[Any] = {
    implicit val actorSystem = cromwellSystem.actorSystem
    implicit val materializer = cromwellSystem.materializer
    implicit val ec = actorSystem.dispatcher
    actorSystem.actorOf(CromwellServerActor.props(cromwellSystem), "cromwell-service")
    Future {
      Await.result(actorSystem.whenTerminated, Duration.Inf)
    }
  }
}

class CromwellServerActor(cromwellSystem: CromwellSystem)(override implicit val materializer: ActorMaterializer)
  extends CromwellRootActor
  with CromwellApiService
  with SwaggerService {
  implicit val actorSystem = context.system
  override implicit val ec = context.dispatcher
  override def actorRefFactory: ActorContext = context

  override val serverMode = true
  override val abortJobsOnTerminate = false

  val webserviceConf = cromwellSystem.conf.getConfig("webservice")
  val interface = webserviceConf.getString("interface")
  val port = webserviceConf.getInt("port")

  val routeUnwrapped: Boolean = cromwellSystem.conf.as[Option[Boolean]]("api.routeUnwrapped").getOrElse(false)

  val allRoutes: Route = routes.wrapped("api", routeUnwrapped) ~ swaggerUiResourceRoute

  Http().bindAndHandle(allRoutes, interface, port) onComplete {
    case Success(_) => actorSystem.log.info("Cromwell service started...")
    case Failure(e) =>
      /*
        TODO:
        If/when CromwellServer behaves like a better async citizen, we may be less paranoid about our async log messages
        not appearing due to the actor system shutdown. For now, synchronously print to the stderr so that the user has
        some idea of why the server failed to start up.
      */
      Console.err.println(s"Binding failed interface $interface port $port")
      e.printStackTrace(Console.err)
      cromwellSystem.shutdownActorSystem()
  }

  /*
    During testing it looked like not explicitly invoking the WMA in order to evaluate all of the lazy actors in
    CromwellRootActor would lead to weird timing issues the first time it was invoked organically
   */
  workflowManagerActor
}

object CromwellServerActor {
  def props(cromwellSystem: CromwellSystem)(implicit materializer: ActorMaterializer): Props = {
    Props(new CromwellServerActor(cromwellSystem)).withDispatcher(EngineDispatcher)
  }
}
