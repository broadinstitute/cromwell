package common.util

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class VersionUtilSpec extends AnyFlatSpec with Matchers {

  behavior of "VersionUtil"

  private val VersionRegex = """([0-9]+)(-[0-9a-f]{7}(-SNAP)?$)?""".r

  it should "getVersion" in {
    val version = VersionUtil.getVersion("cromwell-common")
    version match {
      case VersionRegex(_, _, _) =>
      case "cromwell-common-version.conf-to-be-generated-by-sbt" =>
      case _ => fail(s"Did the versioning scheme change? Got $version")
    }
  }

  it should "getVersion with the default" in {
    val version = VersionUtil.getVersion("made-up-artifact")
    version should be ("made-up-artifact-version.conf-to-be-generated-by-sbt")
  }

  it should "getVersion with a default override" in {
    val version = VersionUtil.getVersion("made-up-artifact", _ => "default override")
    version should be ("default override")
  }

  it should "defaultMessage" in {
    VersionUtil.defaultMessage("some-project") should be("some-project-version.conf-to-be-generated-by-sbt")
  }

  it should "fail sbtDependencyVersion for a made up version" in {
    val expected = intercept[RuntimeException](VersionUtil.sbtDependencyVersion("madeUp")("made-up-project"))
    expected.getMessage should fullyMatch regex
      "Did not parse a version for 'madeUpV' from .*/project/Dependencies.scala " +
        "\\(This occurred after made-up-project-version.conf was not found.\\)"
  }

  it should "pass sbtDependencyVersion check for typesafeConfig" in {
    val version = VersionUtil.sbtDependencyVersion("typesafeConfig")("made-up-typesafe-config-project")
    version should fullyMatch regex """\d+\.\d+.\d+"""
  }

}
