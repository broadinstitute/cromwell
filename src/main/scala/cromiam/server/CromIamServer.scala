package cromiam.server

import akka.Done
import akka.actor.ActorSystem
import akka.event.Logging
import akka.http.scaladsl.server.{HttpApp, Route}
import akka.stream.ActorMaterializer
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.ConfigFactory
import cromiam.server.config.CromIamServerConfig
import cromiam.webservice.{CromIamApiService, SwaggerService}

import scala.concurrent.{ExecutionContext, ExecutionContextExecutor, Future, Promise}


object CromIamServer extends HttpApp with CromIamApiService with SwaggerService {

  final val configuration: CromIamServerConfig = CromIamServerConfig.getFromConfig(ConfigFactory.load()) match {
    case Valid(c) => c
    case Invalid(errors) => throw new Exception("Bad CromIAM configuration:" + errors.toList.mkString("\n", "\n", "\n"))
  }

  def run(): Unit = {
    CromIamServer.startServer(configuration.cromIamConfig.http.interface, configuration.cromIamConfig.http.port, configuration.cromIamConfig.serverSettings)
  }

  override implicit val system: ActorSystem = ActorSystem()
  override implicit lazy val executor: ExecutionContextExecutor = system.dispatcher
  override implicit val materializer: ActorMaterializer = ActorMaterializer()



  override val logger = Logging(system, getClass)

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
