package cromiam.server

import akka.actor.ActorSystem
import akka.event.Logging
import akka.http.scaladsl.Http
import akka.http.scaladsl.server.Route
import akka.stream.Materializer
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.{Config, ConfigFactory}
import common.util.VersionUtil
import cromiam.server.config.{CromIamServerConfig, SwaggerOauthConfig}
import cromiam.server.status.StatusService
import cromiam.webservice.{CromIamApiService, SwaggerService}
import org.broadinstitute.dsde.workbench.util.health.Subsystems.{Cromwell, Sam}

import scala.concurrent.ExecutionContextExecutor

object CromIamServer extends CromIamApiService with SwaggerService {

  final val rootConfig: Config = ConfigFactory.load()

  final val configuration: CromIamServerConfig = CromIamServerConfig.getFromConfig(rootConfig) match {
    case Valid(c) => c
    case Invalid(errors) => throw new Exception("Bad CromIAM configuration:" + errors.toList.mkString("\n", "\n", "\n"))
  }

  override val oauthConfig: SwaggerOauthConfig = configuration.swaggerOauthConfig

  implicit val system: ActorSystem = ActorSystem()
  implicit lazy val executor: ExecutionContextExecutor = system.dispatcher
  implicit val materializer: Materializer = Materializer(system)

  val log = Logging(system, getClass)

  val routes: Route = allRoutes ~ swaggerUiResourceRoute

  override val statusService: StatusService = new StatusService(() =>
    Map(Cromwell -> cromwellClient.subsystemStatus(), Sam -> samClient.subsystemStatus())
  )

  def run(): Unit = {
    log.info(s"Version {}", VersionUtil.getVersion("cromiam"))
    val bindingFuture = Http().newServerAt(
      configuration.cromIamConfig.http.interface,
      configuration.cromIamConfig.http.port)
      .withSettings(configuration.cromIamConfig.serverSettings)
      .bind(routes)

    // Add shutdown hook
    val _ = sys.addShutdownHook {
      log.info("Shutting down the server")
      bindingFuture
        .flatMap(_.unbind())
        .onComplete(_ => system.terminate())
    }
  }
}
