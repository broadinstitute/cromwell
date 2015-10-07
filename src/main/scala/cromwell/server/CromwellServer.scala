package cromwell.server

import akka.io.IO
import akka.pattern.ask
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import cromwell.webservice.{CromwellApiServiceActor, SwaggerService}
import spray.can.Http

import scala.concurrent.duration._
import scala.util.{Failure, Success}

// Note that as per the language specification, this is instantiated lazily and only used when necessary (i.e. server mode)
object CromwellServer extends DefaultWorkflowManagerSystem {
  val conf = ConfigFactory.load()

  val service = actorSystem.actorOf(CromwellApiServiceActor.props(workflowManagerActor, SwaggerService.from(conf)), "cromwell-service")

  implicit val timeout = Timeout(5.seconds)

  val webserviceConf = conf.getConfig("webservice")

  import scala.concurrent.ExecutionContext.Implicits.global
  (IO(Http) ? Http.Bind(service, interface =  webserviceConf.getString("interface"), port = webserviceConf.getInt("port"))).onComplete {
    case Success(Http.CommandFailed(failure)) =>
      actorSystem.log.error("could not bind to port: " + failure.toString)
      actorSystem.shutdown()
    case Failure(t) =>
      actorSystem.log.error(t, "could not bind to port")
      actorSystem.shutdown()
    case _ =>
      actorSystem.log.info("Cromwell service started...")
  }
}
