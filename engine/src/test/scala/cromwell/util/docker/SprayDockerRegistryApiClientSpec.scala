package cromwell.util.docker

import akka.actor.ActorSystem
import com.google.api.client.auth.oauth2.Credential
import com.typesafe.config.ConfigFactory
import cromwell.CromwellSpec.{DockerTest, IntegrationTest}
import cromwell.engine.backend.BackendConfiguration
import cromwell.engine.io.gcs.GoogleConfiguration
import cromwell.util.DockerConfiguration
import cromwell.util.google.{GoogleCredentialFactory, GoogleCredentialFactorySpec}
import org.scalatest.concurrent.{IntegrationPatience, ScalaFutures}
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import spray.http.HttpResponse
import spray.httpx.UnsuccessfulResponseException

import scala.language.postfixOps
import scala.util.Try

class SprayDockerRegistryApiClientSpec extends FlatSpec with Matchers with BeforeAndAfterAll with ScalaFutures
with IntegrationPatience {

  private var actorSystem: ActorSystem = _
  private var client: SprayDockerRegistryApiClient = _

  private val config = BackendConfiguration.DefaultBackendEntry.config
  private val googleConfiguration = GoogleConfiguration.fromConfig(config)
  private val cromwellAuthCredential: Option[Credential] = googleConfiguration flatMap { gc => Try(GoogleCredentialFactory.fromCromwellAuthScheme(gc)) } toOption

  override protected def beforeAll() = {
    actorSystem = ActorSystem("spray-docker-test-actor-system")
    client = new SprayDockerRegistryApiClient()(actorSystem)
  }

  override protected def afterAll() = {
    actorSystem.shutdown()
    actorSystem = null
    client = null
  }

  behavior of "SprayDockerRegistryApiClient"

  it should "resolve docker hub image hashes" taggedAs IntegrationTest in {
    val identifiers = Table(
      "identifier",
      "ubuntu",
      "library/ubuntu",
      "ubuntu:latest",
      "library/ubuntu:latest")

    forAll(identifiers) { identifier =>
      val dockerHash = client.getDockerHash(config, cromwellAuthCredential, identifier).futureValue
      dockerHash.hashType should be("layerBlobs-sha256-md5")
      dockerHash.hashString should have length 32
    }
  }

  it should "resolve docker hub image v1 layer ids" taggedAs IntegrationTest in {
    val identifiers = Table(
      "identifier",
      "ubuntu",
      "library/ubuntu",
      "ubuntu:latest",
      "library/ubuntu:latest")

    forAll(identifiers) { identifier =>
      val parsed = DockerIdentifierParser(config, cromwellAuthCredential).parse(identifier)
      val tagId = parsed.asInstanceOf[DockerTagIdentifier]
      val dockerHash = client.getImageId(tagId).futureValue.dockerHash.get
      dockerHash.hashType should be("layerIds-sha256Part-md5")
      dockerHash.hashString should have length 32
    }
  }

  it should "resolve docker hub and gcr digests" taggedAs IntegrationTest in {
    val identifiers = Table(
      "identifier",
      "ubuntu@sha256:f91f9bab1fe6d0db0bfecc751d127a29d36e85483b1c68e69a246cf1df9b4251",
      "library/ubuntu@sha256:f91f9bab1fe6d0db0bfecc751d127a29d36e85483b1c68e69a246cf1df9b4251",
      // Haven't actually seen a gcr digest yet, so making one up for this test.
      "gcr.io/broad-dsde-dev/fauxbuntu@sha256:f91f9bab1fe6d0db0bfecc751d127a29d36e85483b1c68e69a246cf1df9b4251")

    forAll(identifiers) { identifier =>
      val dockerHash = client.getDockerHash(config, cromwellAuthCredential, identifier).futureValue
      dockerHash.hashType should be("digest-sha256")
      dockerHash.hashString should be("f91f9bab1fe6d0db0bfecc751d127a29d36e85483b1c68e69a246cf1df9b4251")
    }
  }

  it should "resolve docker hub tags using authentication" taggedAs(DockerTest, IntegrationTest) in {
    DockerHubLoginProviderSpec.assumeDockerHubAuthExists()

    val identifiers = Table(
      "identifier",
      "ubuntu",
      "library/ubuntu",
      "ubuntu:latest",
      "library/ubuntu:latest")

    val dockerConf = DockerConfiguration.build(DockerHubLoginProviderSpec.DockerHubConfig)
    val parser = new DockerIdentifierParser(dockerConf, None)

    forAll(identifiers) { identifier =>
      val parsed = parser.parse(identifier)
      val dockerHash = client.getDockerHashable(parsed).futureValue.dockerHash.get
      dockerHash.hashType should be("layerBlobs-sha256-md5")
      dockerHash.hashString should have length 32
    }
  }

  it should "not resolve docker hub tags when using bad authentication" taggedAs IntegrationTest in {
    val identifiers = Table(
      "identifier",
      "ubuntu",
      "library/ubuntu",
      "ubuntu:latest",
      "library/ubuntu:latest")

    val badDockerConfig = ConfigFactory.parseString(
      s"""
         |docker {
         |  dockerAccount = "fakeAccount"
         |  dockerToken = "YmFkdXNlcjpiYWRwYXNz" // baduser:badpass
         |}
     """.stripMargin)

    val dockerConf = DockerConfiguration.build(badDockerConfig)
    val parser = new DockerIdentifierParser(dockerConf, None)

    forAll(identifiers) { identifier =>
      val parsed = parser.parse(identifier)
      val exception = client.getDockerHashable(parsed).failed.futureValue
      exception shouldBe a[UnsuccessfulResponseException]
      exception.getMessage should startWith("Status: 401 Unauthorized\nBody:")
    }
  }

  it should "fail to resolve gcr image tags when google authentication is not setup" taggedAs IntegrationTest in {
    val identifiers = Table(
      "identifier",
      "gcr.io/broad-dsde-dev/ubuntu",
      "us.gcr.io/broad-dsde-dev/cromwell:dev")

    forAll(identifiers) { identifier =>
      val exception = client.getDockerHash(config, cromwellAuthCredential, identifier).failed.futureValue
      exception shouldBe an[UnsuccessfulResponseException]
      exception.getMessage should startWith("Status: 404 Not Found\nBody:")
    }
  }

  it should "fail to resolve gcr image tags when unauthenticated" taggedAs IntegrationTest in {
    val gcrRegistry = DockerRegistry("gcr.io", NoLoginProvider)
    val gcrUSRegistry = DockerRegistry("us.gcr.io", NoLoginProvider)

    val identifiers = Table(
      "identifier",
      DockerTagIdentifier("broad-dsde-dev/ubuntu", "latest", gcrRegistry),
      DockerTagIdentifier("broad-dsde-dev/cromwell", "dev", gcrUSRegistry))

    forAll(identifiers) { identifier =>
      val exception = client.getDockerHashable(identifier).failed.futureValue
      exception shouldBe an[UnsuccessfulResponseException]
      /*
        Something on the GCR server changed, or perhaps our test environment now supplies a default credential? Either
        way, we're now getting 403 Forbidden errors instead of 404 Not Found. Just in case this a temporary change,
        check for either. Alternatively, we could also just check for any 4xx, or not test for the message at all.
       */
      exception.getMessage should (
        be(s"Status: 403 Forbidden\nBody: Unable to access the repository: ${identifier.name}; " +
          "please verify that it exists and you have permission to access it (no valid credential was supplied).") or
        be("Status: 404 Not Found\nBody: Not found."))
    }
  }

  it should "resolve gcr image tags when google authentication is setup" taggedAs IntegrationTest in {
    GoogleCredentialFactorySpec.assumeAccountConfigExists()

    val identifiers = Table(
      "identifier",
      "gcr.io/broad-dsde-dev/ubuntu",
      "us.gcr.io/broad-dsde-dev/cromwell:dev")

    val parser = DockerIdentifierParser(config, Option(GoogleCredentialFactory.fromCromwellAuthScheme(GoogleCredentialFactorySpec.GoogleAccountConfig)))

    forAll(identifiers) { identifier =>
      val parsed = parser.parse(identifier)
      val dockerHash = client.getDockerHashable(parsed).futureValue.dockerHash.get
      dockerHash.hashType should be("imageId-sha256")
      dockerHash.hashString should have length 64
    }
  }

  it should "throw an error when a response doesn't contain docker v1 token header" in {
    val exception = intercept[DockerHeaderNotFoundException](client.toDockerV1Token(HttpResponse()))
    exception.getMessage should be("Response did not contain header X-Docker-Token")
  }

  it should "throw an error when a response doesn't contain docker v2 token header" in {
    val exception = intercept[DockerHeaderNotFoundException](client.toDockerV2TokenRequest(HttpResponse()))
    exception.getMessage should be("Response did not contain header WWW-Authenticate")
  }
}

object NoLoginProvider extends DockerLoginProvider {
  override def dockerLogin = None
}
