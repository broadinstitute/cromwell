package server

import akka.Done
import akka.actor.ActorSystem
import akka.event.Logging
import akka.http.scaladsl.server.{HttpApp, Route}
import akka.http.scaladsl.settings.ServerSettings
import akka.stream.ActorMaterializer
import com.typesafe.config.{Config, ConfigFactory}
import webservice.{CromIamApiService, SwaggerService}

import scala.concurrent.{ExecutionContext, ExecutionContextExecutor, Future, Promise}


object CromIamServer extends HttpApp with CromIamApiService with SwaggerService {

  def run(): Unit = {
    CromIamServer.startServer(CromIamServer.config.getString("http.interface"), CromIamServer.config.getInt("http.port"), ServerSettings(CromIamServer.config))
  }

  override implicit val system: ActorSystem = ActorSystem()
  override implicit lazy val executor: ExecutionContextExecutor = system.dispatcher
  override implicit val materializer: ActorMaterializer = ActorMaterializer()

  private val config: Config = ConfigFactory.load()
  override final val cromwellInterface: String = config.getString("cromwell.interface")
  override final val cromwellPort: Int = config.getInt("cromwell.port")


  val logger = Logging(system, getClass)
  override val route: Route = allRoutes ~ swaggerUiResourceRoute

  // Override default shutdownsignal which was just "hit return/enter"
  override def waitForShutdownSignal(actorSystem: ActorSystem)(implicit executionContext: ExecutionContext): Future[Done] = {
    val promise = Promise[Done]()
    sys.addShutdownHook {
      // we can add anything we want the server to do when someone shutdowns the server (Ctrl-c)
      logger.info("Shutting down the server")
      promise.success(Done)
    }
    promise.future
  }
}
