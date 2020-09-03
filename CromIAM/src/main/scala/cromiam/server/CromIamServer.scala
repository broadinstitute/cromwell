package cromiam.server

import akka.Done
import akka.actor.ActorSystem
import akka.event.Logging
import akka.http.scaladsl.server.{HttpApp, Route}
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.{Config, ConfigFactory}
import common.util.VersionUtil
import cromiam.server.config.{CromIamServerConfig, SwaggerOauthConfig}
import cromiam.server.status.StatusService
import cromiam.webservice.{CromIamApiService, SwaggerService}
import org.broadinstitute.dsde.workbench.util.health.Subsystems.{Cromwell, Sam}

import scala.concurrent.{ExecutionContext, ExecutionContextExecutor, Future, Promise}


object CromIamServer extends HttpApp with CromIamApiService with SwaggerService {

  final val rootConfig: Config = ConfigFactory.load()

  final val configuration: CromIamServerConfig = CromIamServerConfig.getFromConfig(rootConfig) match {
    case Valid(c) => c
    case Invalid(errors) => throw new Exception("Bad CromIAM configuration:" + errors.toList.mkString("\n", "\n", "\n"))
  }

  override val oauthConfig: SwaggerOauthConfig = configuration.swaggerOauthConfig

  def run(): Unit = {
    log.info(s"Version {}", VersionUtil.getVersion("cromiam"))
    /*
    TODO: We're not passing in our actor system so a new one is getting created by akka.http.scaladsl.server.HttpApp.
    Things seem stable at the moment so not touching it. Should investigate if we should be sharing the ActorSystem.
    If there is a reason then leave a comment why there should be two actor systems.
    https://github.com/broadinstitute/cromwell/issues/3851
     */
    CromIamServer.startServer(configuration.cromIamConfig.http.interface, configuration.cromIamConfig.http.port, configuration.cromIamConfig.serverSettings)
  }

  override implicit val system: ActorSystem = ActorSystem()
  override implicit lazy val executor: ExecutionContextExecutor = system.dispatcher

  override val log = Logging(system, getClass)

  override val routes: Route = allRoutes ~ swaggerUiResourceRoute

  override val statusService: StatusService = new StatusService(() => Map(Cromwell -> cromwellClient.subsystemStatus, Sam -> samClient.subsystemStatus))

  // Override default shutdownsignal which was just "hit return/enter"
  override def waitForShutdownSignal(actorSystem: ActorSystem)(implicit executionContext: ExecutionContext): Future[Done] = {
    val promise = Promise[Done]()
    sys.addShutdownHook {
      // we can add anything we want the server to do when someone shutdowns the server (Ctrl-c)
      log.info("Shutting down the server")
      promise.success(Done)
    }
    promise.future
  }
}
