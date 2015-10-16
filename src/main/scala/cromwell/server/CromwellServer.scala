package cromwell.server

import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import cromwell.webservice.{CromwellApiServiceActor, SwaggerService}
import lenthall.spray.SprayCanHttpService._

import scala.concurrent.duration._

// Note that as per the language specification, this is instantiated lazily and only used when necessary (i.e. server mode)
object CromwellServer extends DefaultWorkflowManagerSystem {
  val conf = ConfigFactory.load()

  val service = actorSystem.actorOf(CromwellApiServiceActor.props(workflowManagerActor, SwaggerService.from(conf)), "cromwell-service")

  implicit val timeout = Timeout(5.seconds)

  val webserviceConf = conf.getConfig("webservice")

  import scala.concurrent.ExecutionContext.Implicits.global

  service.bindOrShutdown(
    interface = webserviceConf.getString("interface"),
    port = webserviceConf.getInt("port")
  ) onSuccess {
    case _ =>
      actorSystem.log.info("Cromwell service started...")
  }
}
