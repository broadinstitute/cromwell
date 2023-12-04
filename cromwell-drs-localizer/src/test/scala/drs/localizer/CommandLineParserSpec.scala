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
  private val azureIdentityClientId = "itme@azure.com"
  private val manifestPath = "/my/manifest.txt"

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
    args.azureIdentityClientId shouldBe empty
    args.manifestPath shouldBe empty
  }

  it should "successfully parse with three arguments" in {
    val args = parser.parse(Array(drsObject, containerPath, requesterPaysProject), CommandLineArguments()).get
    args.drsObject.get shouldBe drsObject
    args.containerPath.get shouldBe containerPath
    args.accessTokenStrategy.get shouldBe AccessTokenStrategy.Google
    args.googleRequesterPaysProject.get shouldBe requesterPaysProject
    args.azureIdentityClientId shouldBe empty
    args.manifestPath shouldBe empty
  }

  it should "successfully parse requester pays project" in {
    val args = parser.parse(Array(drsObject, containerPath, "-r", requesterPaysProject), CommandLineArguments()).get

    args.drsObject.get shouldBe drsObject
    args.containerPath.get shouldBe containerPath
    args.accessTokenStrategy.get shouldBe AccessTokenStrategy.Google
    args.googleRequesterPaysProject.get shouldBe requesterPaysProject
    args.azureIdentityClientId shouldBe empty
    args.manifestPath shouldBe empty
  }

  it should "successfully parse with three arguments and requester pays project" in {
    val args = parser
      .parse(Array(drsObject, containerPath, requesterPaysProject, "-r", requesterPaysProject), CommandLineArguments())
      .get

    args.drsObject.get shouldBe drsObject
    args.containerPath.get shouldBe containerPath
    args.accessTokenStrategy.get shouldBe AccessTokenStrategy.Google
    args.googleRequesterPaysProject.get shouldBe requesterPaysProject
    args.azureIdentityClientId shouldBe empty
    args.manifestPath shouldBe empty
  }

  it should "fail if requester pays argument and flag specify different projects" in {
    parser.parse(Array(drsObject, containerPath, requesterPaysProject, "-r", "boom!"),
                 CommandLineArguments()
    ) shouldBe None
  }

  it should "successfully parse args with a manifest file" in {
    val args = parser.parse(Array("-m", manifestPath), CommandLineArguments()).get

    args.drsObject shouldBe empty
    args.containerPath shouldBe empty
    args.accessTokenStrategy.get shouldBe AccessTokenStrategy.Google
    args.googleRequesterPaysProject shouldBe empty
    args.azureIdentityClientId shouldBe empty
    args.manifestPath.get shouldBe manifestPath
  }

  it should "fail to parse with a manifest file and one single-file arg" in {
    val args = parser.parse(Array(drsObject, "--manifest-path", manifestPath), CommandLineArguments())
    args shouldBe None
  }

  it should "fail to parse with a manifest file and two single-file args" in {
    val args = parser.parse(Array(drsObject, containerPath, "--manifest-path", manifestPath), CommandLineArguments())
    args shouldBe None
  }

  it should "successfully parse an explicit Google access token strategy invocation" in {
    val args = parser
      .parse(Array(
               "--access-token-strategy",
               "google",
               drsObject,
               containerPath,
               "--requester-pays-project",
               requesterPaysProject
             ),
             CommandLineArguments()
      )
      .get

    args.drsObject.get shouldBe drsObject
    args.containerPath.get shouldBe containerPath
    args.accessTokenStrategy.get shouldBe AccessTokenStrategy.Google
    args.googleRequesterPaysProject.get shouldBe requesterPaysProject
    args.azureIdentityClientId shouldBe empty
    args.manifestPath shouldBe empty
  }

  it should "fail to parse an Azure invocation that specifies requester pays" in {
    val args = parser.parse(
      Array("--access-token-strategy",
            AccessTokenStrategy.Azure,
            drsObject,
            containerPath,
            "--requester-pays-project",
            requesterPaysProject
      ),
      CommandLineArguments()
    )

    args shouldBe None
  }

  it should "successfully parse an Azure invocation" in {
    val args = parser
      .parse(Array("--access-token-strategy", AccessTokenStrategy.Azure, drsObject, containerPath),
             CommandLineArguments()
      )
      .get

    args.drsObject.get shouldBe drsObject
    args.containerPath.get shouldBe containerPath
    args.accessTokenStrategy.get shouldBe AccessTokenStrategy.Azure
    args.googleRequesterPaysProject shouldBe empty
    args.azureIdentityClientId shouldBe empty
    args.manifestPath shouldBe empty
  }

  it should "successfully parse an Azure invocation with identity" in {
    val args = parser
      .parse(
        Array("--access-token-strategy",
              AccessTokenStrategy.Azure,
              "--identity-client-id",
              azureIdentityClientId,
              drsObject,
              containerPath
        ),
        CommandLineArguments()
      )
      .get

    args.drsObject.get shouldBe drsObject
    args.containerPath.get shouldBe containerPath
    args.accessTokenStrategy.get shouldBe AccessTokenStrategy.Azure
    args.googleRequesterPaysProject shouldBe empty
    args.azureIdentityClientId.get shouldBe azureIdentityClientId
    args.manifestPath shouldBe empty
  }

  it should "fail to parse with an unrecognized access token strategy" in {
    val args =
      parser.parse(Array("--access-token-strategy", "nebulous", drsObject, containerPath), CommandLineArguments())
    args shouldBe None
  }
}
