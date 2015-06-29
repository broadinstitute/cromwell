package cromwell.server

import java.io.File

import akka.io.IO
import akka.pattern.ask
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import com.wordnik.swagger.model.ApiInfo
import cromwell.engine.db.DataAccess
import cromwell.engine.db.slick.DataAccessController
import cromwell.webservice.{CromwellApiService, CromwellApiServiceActor, SwaggerService}
import spray.can.Http

import scala.concurrent.duration._
import scala.reflect.runtime.universe._
import scala.util.{Failure, Success}

// Note that as per the language specification, this is instantiated lazily and only used when necessary (i.e. server mode)
object CromwellServer extends DefaultWorkflowManagerSystem {
  val conf = ConfigFactory.parseFile(new File("/etc/cromwell.conf"))
  private lazy val realDataAccess = DataAccessController

  val swaggerConfig = conf.getConfig("swagger")
  val swaggerService = new SwaggerService(
    swaggerConfig.getString("apiVersion"),
    swaggerConfig.getString("baseUrl"),
    swaggerConfig.getString("apiDocs"),
    swaggerConfig.getString("swaggerVersion"),
    Vector(typeOf[CromwellApiService]),
    Option(new ApiInfo(
      swaggerConfig.getString("info"),
      swaggerConfig.getString("description"),
      swaggerConfig.getString("termsOfServiceUrl"),
      swaggerConfig.getString("contact"),
      swaggerConfig.getString("license"),
      swaggerConfig.getString("licenseUrl"))
    ))

  val service = actorSystem.actorOf(CromwellApiServiceActor.props(workflowManagerActor, swaggerService), "cromwell-service")

  implicit val timeout = Timeout(5.seconds)

  import scala.concurrent.ExecutionContext.Implicits.global
  (IO(Http) ? Http.Bind(service, interface =  conf.getString("webservice.interface"), port = conf.getInt("webservice.port"))).onComplete {
    case Success(Http.CommandFailed(failure)) =>
      actorSystem.log.error("could not bind to port: " + failure.toString)
      actorSystem.shutdown()
    case Failure(t) =>
      actorSystem.log.error(t, "could not bind to port")
      actorSystem.shutdown()
    case _ =>
      actorSystem.log.info("Cromwell service started...")
  }

  override def dataAccess: DataAccess = realDataAccess
}

