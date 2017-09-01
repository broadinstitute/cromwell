package cromwell.services.healthmonitor.impl

import cromwell.core.docker.{DockerCliClient, DockerCliKey}
import cromwell.services.EngineServicesStore
import cromwell.services.healthmonitor.HealthMonitorServiceActor.{MonitoredSubsystem, OkStatus, SubsystemStatus}

import scala.concurrent.{ExecutionContext, Future}

trait CommonMonitoredSubsystems {
  implicit val ec: ExecutionContext

  lazy val DockerHub = MonitoredSubsystem("DockerHub", checkDockerhub _)
  lazy val EngineDb = MonitoredSubsystem("Engine Database", checkEngineDb _)

  /**
    * Demonstrates connectivity to DockerHub by periodically pulling a small image. Remove the image afterwards in
    * order to ensure that subsequent attempts aren't showing a false positive
    */
  private def checkDockerhub(): Future[SubsystemStatus] = {
    println("DOCKER")
 //   Future {
      // Using a tiny docker image the user is very unlikely to be actively using
      println("YO")
      val alpine = DockerCliKey("alpine", "latest") // FIXME: make a clone of alpine in our dockerhub and use that
      println("HO")

    val res = DockerCliClient.pull(alpine) flatMap { _ => DockerCliClient.rmi(alpine) }

//      val res = for {
//        p <- DockerCliClient.pull(alpine)
//        r <- DockerCliClient.rmi(alpine)
//      } yield r

      println("FOO: " + res)

//      res.get
 //   } map { _ => OkStatus }

    Future.fromTry(res) map { _ => OkStatus }
  }

  /**
    * Demonstrates connectivity to the engine database by periodically making a small query
    */
  private def checkEngineDb(): Future[SubsystemStatus] = {
    println("DB")
    EngineServicesStore.engineDatabaseInterface.queryDockerHashStoreEntries("DOESNOTEXIST") map { _ => OkStatus }
  }
}
