package cromwell.server

import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import cromwell.webservice.CromwellApiServiceActor
import lenthall.spray.SprayCanHttpService._

import scala.concurrent.duration._

// Note that as per the language specification, this is instantiated lazily and only used when necessary (i.e. server mode)
object CromwellServer extends WorkflowManagerSystem {
  implicit val timeout = Timeout(5.seconds)
  import scala.concurrent.ExecutionContext.Implicits.global

  val conf = ConfigFactory.load()
  val service = actorSystem.actorOf(CromwellApiServiceActor.props(workflowManagerActor, conf), "cromwell-service")
  val webserviceConf = conf.getConfig("webservice")
  service.bindOrShutdown(interface = webserviceConf.getString("interface"), port = webserviceConf.getInt("port")) onSuccess {
    case _ => actorSystem.log.info("Cromwell service started...")
  }
}
