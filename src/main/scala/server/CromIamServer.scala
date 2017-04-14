package server

import akka.Done
import akka.actor.ActorSystem
import akka.event.Logging
import akka.http.scaladsl.server.{HttpApp, Route}
import akka.http.scaladsl.settings.ServerSettings
import akka.stream.ActorMaterializer
import com.typesafe.config.ConfigFactory
import service.CromIamApiService

import scala.concurrent.{ExecutionContext, Future, Promise}


object CromIamServer extends HttpApp with CromIamApiService {

  def run() = {
    CromIamServer.startServer(CromIamServer.config.getString("http.interface"), CromIamServer.config.getInt("http.port"), ServerSettings(CromIamServer.config))
  }

  implicit val system = ActorSystem()
  implicit val executor = system.dispatcher
  implicit val materializer = ActorMaterializer()

  val config = ConfigFactory.load()
  val logger = Logging(system, getClass)
  val route: Route = allRoutes

  // Override default shutdownsignal which was just "hit return/enter"
  override def waitForShutdownSignal(actorSystem: ActorSystem)(implicit executionContext: ExecutionContext): Future[Done] = {
    val promise = Promise[Done]()
    sys.addShutdownHook {
      promise.success(Done)
    }
    promise.future
  }
}
