package drs.localizer

import common.assertion.CromwellTimeoutSpec
import drs.localizer.CommandLineParser.AccessTokenStrategy
import org.scalatest.BeforeAndAfter
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class CommandLineParserSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with BeforeAndAfter {

  var parser: scopt.OptionParser[CommandLineArguments] = _

  private val drsObject = "drs://provider/object"
  private val containerPath = "/cromwell_root/my.bam"
  private val requesterPaysProject = "i_heart_egress"
  private val azureVaultName = "Kwikset"
  private val azureSecretName = "shhh"
  private val azureIdentityClientId = "itme@azure.com"

  behavior of "DRS Localizer command line parser"

  before {
    parser = DrsLocalizerMain.buildParser()
  }

  it should "fail to parse with no arguments" in {
    parser.parse(Array.empty[String], CommandLineArguments()) shouldBe None
  }

  it should "fail to parse with only one argument" in {
    parser.parse(Array("some arg"), CommandLineArguments()) shouldBe None
  }

  it should "successfully parse with two arguments" in {
    val args = parser.parse(Array(drsObject, containerPath), CommandLineArguments()).get

    args.drsObject.get shouldBe drsObject
    args.containerPath.get shouldBe containerPath
    args.accessTokenStrategy.get shouldBe AccessTokenStrategy.Google
    args.googleRequesterPaysProject shouldBe empty
    args.azureVaultName shouldBe empty
    args.azureSecretName shouldBe empty
    args.azureIdentityClientId shouldBe empty
  }

  it should "successfully parse with three arguments" in {
    val args = parser.parse(Array(drsObject, containerPath, requesterPaysProject), CommandLineArguments()).get

    args.drsObject.get shouldBe drsObject
    args.containerPath.get shouldBe containerPath
    args.accessTokenStrategy.get shouldBe AccessTokenStrategy.Google
    args.googleRequesterPaysProject.get shouldBe requesterPaysProject
    args.azureVaultName shouldBe empty
    args.azureSecretName shouldBe empty
    args.azureIdentityClientId shouldBe empty
  }

  it should "successfully parse an explicit Google access token strategy invocation" in {
    val args = parser.parse(Array("--access-token-strategy", "google", drsObject, containerPath, requesterPaysProject), CommandLineArguments()).get

    args.drsObject.get shouldBe drsObject
    args.containerPath.get shouldBe containerPath
    args.accessTokenStrategy.get shouldBe AccessTokenStrategy.Google
    args.googleRequesterPaysProject.get shouldBe requesterPaysProject
    args.azureVaultName shouldBe empty
    args.azureSecretName shouldBe empty
    args.azureIdentityClientId shouldBe empty
  }

  it should "fail to parse an Azure invocation missing vault name and secret name" in {
    val args = parser.parse(Array(
      "--access-token-strategy", AccessTokenStrategy.Azure,
      drsObject, containerPath), CommandLineArguments())

    args shouldBe None
  }

  it should "fail to parse an Azure invocation missing vault name" in {
    val args = parser.parse(Array(
      "--access-token-strategy", AccessTokenStrategy.Azure,
      "--secret-name", azureSecretName,
      drsObject, containerPath), CommandLineArguments())

    args shouldBe None
  }

  it should "fail to parse an Azure invocation missing secret name" in {
    val args = parser.parse(Array(
      "--access-token-strategy", AccessTokenStrategy.Azure,
      "--vault-name", azureVaultName,
      drsObject, containerPath), CommandLineArguments())

    args shouldBe None
  }

  it should "fail to parse an Azure invocation that specifies requester pays" in {
    val args = parser.parse(Array(
      "--access-token-strategy", AccessTokenStrategy.Azure,
      "--secret-name", azureSecretName,
      "--vault-name", azureVaultName,
      drsObject, containerPath, requesterPaysProject), CommandLineArguments())

    args shouldBe None
  }

  it should "successfully parse an Azure invocation" in {
    val args = parser.parse(Array(
      "--access-token-strategy", AccessTokenStrategy.Azure,
      "--secret-name", azureSecretName,
      "--vault-name", azureVaultName,
      drsObject, containerPath), CommandLineArguments()).get

    args.drsObject.get shouldBe drsObject
    args.containerPath.get shouldBe containerPath
    args.accessTokenStrategy.get shouldBe AccessTokenStrategy.Azure
    args.googleRequesterPaysProject shouldBe empty
    args.azureVaultName.get shouldBe azureVaultName
    args.azureSecretName.get shouldBe azureSecretName
    args.azureIdentityClientId shouldBe empty
  }

  it should "successfully parse an Azure invocation with all the trimmings" in {
    val args = parser.parse(Array(
      "--access-token-strategy", AccessTokenStrategy.Azure,
      "--vault-name", azureVaultName,
      "--secret-name", azureSecretName,
      "--identity-client-id", azureIdentityClientId,
      drsObject, containerPath), CommandLineArguments()).get

    args.drsObject.get shouldBe drsObject
    args.containerPath.get shouldBe containerPath
    args.accessTokenStrategy.get shouldBe AccessTokenStrategy.Azure
    args.googleRequesterPaysProject shouldBe empty
    args.azureVaultName.get shouldBe azureVaultName
    args.azureSecretName.get shouldBe azureSecretName
    args.azureIdentityClientId.get shouldBe azureIdentityClientId
  }

  it should "fail to parse with an unrecognized access token strategy" in {
    val args = parser.parse(Array("--access-token-strategy", "nebulous", drsObject, containerPath), CommandLineArguments())
    args shouldBe None
  }
}
