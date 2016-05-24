package cromwell.util.docker

import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.CromwellSpec.{DockerTest, IntegrationTest}
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{Assertions, FlatSpec, Matchers}

import scala.util.Try
import scala.util.parsing.json.{JSON, JSONObject}

class DockerHubLoginProviderSpec extends FlatSpec with Matchers {
  behavior of "DockerHubLoginProvider"

  it should "decode a docker token into username and password" in {
    val tokens = Table(
      ("token", "username", "password"),
      ("aGVsbG86d29ybGQ=", "hello", "world"),
      ("aGVsbG86Ondvcmxk", "hello", ":world"),
      ("aGVsbG86", "hello", ""),
      ("aGVsbG8=", "hello", ""),
      ("Ondvcmxk", "", "world"),
      ("", "", "")
    )

    forAll(tokens) { (token, username, password) =>
      val loginProvider = new DockerHubLoginProvider(Option(token))
      val loginOption = loginProvider.dockerLogin
      loginOption shouldNot be(empty)
      val login = loginOption.get
      login should be(DockerLogin(username, password))
    }
  }

  it should "not provide a username and password for an unspecified token" in {
    val loginProvider = new DockerHubLoginProvider(None)
    val loginOption = loginProvider.dockerLogin
    loginOption should be(empty)
  }

  it should "not parse an invalid docker token" in {
    val exception = intercept[IllegalArgumentException](new DockerHubLoginProvider(Option("bad_base64")))
    exception shouldBe an[IllegalArgumentException]
    exception.getMessage should be("Illegal base64 character 5f")
  }

  it should "find the docker hub token's username and password" taggedAs (DockerTest, IntegrationTest) in {
    DockerHubLoginProviderSpec.assumeDockerHubAuthExists()

    val user = new DockerHubLoginProvider(DockerHubLoginProviderSpec.DockerHubConfig).dockerLogin.get.username
    user shouldNot be(empty)
  }
}

object DockerHubLoginProviderSpec extends {
  val DockerConfigFile = home / ".docker" / "config.json"

  lazy val DockerHubAuthTry = Try {
    val jsonString = DockerConfigFile.contentAsString
    val jsonObject = JSON.parseRaw(jsonString).get.asInstanceOf[JSONObject]
    val jsonAuths = jsonObject.obj("auths").asInstanceOf[JSONObject]
    val dockerAuthObject = jsonAuths.obj("https://index.docker.io/v1/").asInstanceOf[JSONObject]
    val dockerAccount = dockerAuthObject.obj("email").asInstanceOf[String]
    val dockerToken = dockerAuthObject.obj("auth").asInstanceOf[String]
    (dockerAccount, dockerToken)
  }

  lazy val DockerHubAuth = DockerHubAuthTry.get

  lazy val DockerHubAuthExists = DockerHubAuthTry.isSuccess

  lazy val DockerHubConfig = ConfigFactory.parseString(
      s"""
         |docker {
         |  dockerAccount = "${DockerHubAuth._1}"
         |  dockerToken = "${DockerHubAuth._2}"
         |}
     """.stripMargin)

  import Assertions._

  def assumeDockerHubAuthExists() = assume(DockerHubAuthExists,
    s"""
       |Unable to get docker hub auth from $DockerConfigFile
       |${DockerHubAuthTry.failed.toOption.fold("")(_.toString)}
       |Did you run `docker login`?
       |""".stripMargin)
}
