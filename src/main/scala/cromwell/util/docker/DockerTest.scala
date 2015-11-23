package cromwell.util.docker

import akka.actor.ActorSystem

import scala.concurrent.Await
import scala.concurrent.duration.Duration

object DockerTest {
  implicit val actorSystem = ActorSystem("docker-test-actor-system")

  implicit val executionContext = actorSystem.dispatcher

  def run() = {
    val tests = Seq(
      "ubuntu@sha256:f91f9bab1fe6d0db0bfecc751d127a29d36e85483b1c68e69a246cf1df9b4251",
      "broadinstitute/scala-baseimage",
      "us.gcr.io/broad-dsde-dev/cromwell:dev",
      "gcr.io/broad-dsde-dev/ubuntu",
      "ubuntu")

    val client: DockerRegistryApiClient = new SprayDockerRegistryApiClient()

    for (identifier <- tests) {
      val parsed = DockerIdentifierParser.parse(identifier)
      val futureHashable = client.getDockerHashable(parsed) map { hashable =>
        println(
          s"""SUCCESS
              |parsed:   $identifier
              |hashable: $hashable
              |""".stripMargin)
      } recover {
        case error =>
          println(
            s"""FAILED
                |parsed:   $identifier
                |error:    $error
                |""".stripMargin)
      }
      Await.result(futureHashable, Duration.Inf)
    }
  }
}
