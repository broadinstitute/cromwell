package cromwell.util.docker

import akka.actor.ActorSystem

import scala.concurrent.Await
import scala.concurrent.duration.Duration

object DockerTest {
  implicit val actorSystem = ActorSystem("docker-test-actor-system")

  implicit val executionContext = actorSystem.dispatcher

  // TODO: Need to wire this to google's OAuth, like the JES client is currently doing.
  // TODO: For now, generate your own token with "gcloud auth print-access-token".
  val gcrToken = "ya29.0000_0000-1111_1111-2222_2222-3333_3333-4444_4444-5555_5555-6666_6666-7"
  val gcrLogin = Option(DockerLogin("_token", gcrToken))

  def run() = {
    val tests = Map(
      "ubuntu@sha256:f91f9bab1fe6d0db0bfecc751d127a29d36e85483b1c68e69a246cf1df9b4251" -> None,
      "broadinstitute/scala-baseimage" -> None,
      "us.gcr.io/broad-dsde-dev/cromwell:dev" -> gcrLogin,
      "gcr.io/broad-dsde-dev/ubuntu" -> gcrLogin)

    val client: DockerRegistryApiClient = new SprayDockerRegistryApiClient()

    for ((identifier, login) <- tests) {
      val parsed = DockerIdentifierParser.parse(identifier)
      val futureHashable = client.getDockerHashable(parsed, login) map { hashable =>
        println(
          s"""SUCCESS
              |parsed:   $parsed
              |hashable: $hashable
              |""".stripMargin)
      } recover {
        case error =>
          println(
            s"""FAILED
                |parsed:   $parsed
                |error:    $error
                |""".stripMargin)
      }
      Await.result(futureHashable, Duration.Inf)
    }
  }
}
