package cromwell

import java.io.File
import java.nio.file.{Files, Paths}

import akka.actor.{Props, ActorSystem}
import akka.io.IO
import akka.pattern.ask
import akka.util.Timeout
import com.typesafe.config.{Config, ConfigFactory}
import com.wordnik.swagger.model.ApiInfo
import cromwell.engine.WorkflowManagerActor
import cromwell.webservice._
import spray.can.Http

import scala.concurrent.duration._
import scala.reflect.runtime.universe._
import scala.util.{Failure, Success}


object Boot extends App {

  private def startup(): Unit = {
    val conf = ConfigFactory.parseFile(new File("/etc/cromwell.conf"))

    // we need an ActorSystem to host our application in
    implicit val system = ActorSystem("cromwell")

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

    val workflowManagerActorRef = system.actorOf(WorkflowManagerActor.props)

    val service = system.actorOf(CromwellApiServiceActor.props(workflowManagerActorRef, swaggerService), "cromwell-service")

    implicit val timeout = Timeout(5.seconds)

    import scala.concurrent.ExecutionContext.Implicits.global
    (IO(Http) ? Http.Bind(service, interface =  conf.getString("webservice.interface"), port = conf.getInt("webservice.port"))).onComplete {
      case Success(Http.CommandFailed(failure)) =>
        system.log.error("could not bind to port: " + failure.toString)
        system.shutdown()
      case Failure(t) =>
        system.log.error(t, "could not bind to port")
        system.shutdown()
      case _ =>
        system.log.info("Cromwell service started...")
    }
  }

  startup()
}
