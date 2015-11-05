package cromwell.util.docker

import cromwell.util.DockerConfiguration
import cromwell.util.google.{GoogleCredentialFactory, GoogleCredentialFactorySpec}
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}

class DockerIdentifierParserSpec extends FlatSpec with Matchers {
  behavior of "DockerIdentifierParser"

  it should "parse docker tagged identifiers" in {
    val parser = DockerIdentifierParser.Default

    val identifiers = Table(
      ("identifier", "name", "tag", "namespace"),
      ("team/imageName", "team/imageName", "latest", "docker.io"),
      ("team/imageName:tag", "team/imageName", "tag", "docker.io"),
      ("imageName", "library/imageName", "latest", "docker.io"),
      ("imageName:", "library/imageName:", "latest", "docker.io"),
      ("imageName:tag", "library/imageName", "tag", "docker.io"),

      // TODO: We should be able to handle other registries _correctly_.
      // TODO: For now, assume the host is actually a user/team.
      ("quay.io/namespace/repository", "quay.io/namespace/repository", "latest", "docker.io"),
      ("quay.io/namespace/repository:tag", "quay.io/namespace/repository", "tag", "docker.io"))

    forAll(identifiers) { (identifier, name, tag, namespace) =>
      val parsed = parser.parse(identifier)
      parsed shouldBe a[DockerTagIdentifier]
      val tagged = parsed.asInstanceOf[DockerTagIdentifier]
      tagged.name should be(name)
      tagged.tag should be(tag)
      tagged.registry.namespace should be(namespace)
    }
  }

  it should "parse gcr.io tagged identifiers" in {
    GoogleCredentialFactorySpec.assumeAccountConfigExists()

    val googleCredentials = new GoogleCredentialFactory {
      override val GoogleConf = GoogleCredentialFactorySpec.GoogleAccountConfig
    }.fromCromwellAuthScheme

    val dockerConf = DockerConfiguration.build(GoogleCredentialFactorySpec.AccountConfig)

    val parser = new DockerIdentifierParser(dockerConf, Option(googleCredentials))

    val identifiers = Table(
      ("identifier", "name", "tag", "namespace"),
      ("gcr.io/google-project/imageName", "google-project/imageName", "latest", "gcr.io"),
      ("gcr.io/google-project/imageName:tag", "google-project/imageName", "tag", "gcr.io"),
      ("us.gcr.io/google-project/imageName", "google-project/imageName", "latest", "us.gcr.io"),
      ("eu.gcr.io/google-project/imageName", "google-project/imageName", "latest", "eu.gcr.io"),
      ("asia.gcr.io/google-project/imageName", "google-project/imageName", "latest", "asia.gcr.io"),
      ("b.gcr.io/google-bucket/imageName", "google-bucket/imageName", "latest", "b.gcr.io"))

    forAll(identifiers) { (identifier, name, tag, namespace) =>
      val parsed = parser.parse(identifier)
      parsed shouldBe a[DockerTagIdentifier]
      val tagged = parsed.asInstanceOf[DockerTagIdentifier]
      tagged.name should be(name)
      tagged.tag should be(tag)
      tagged.registry.namespace should be(namespace)
    }
  }
}
