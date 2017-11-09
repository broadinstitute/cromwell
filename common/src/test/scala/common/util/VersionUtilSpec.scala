package common.util

import org.scalatest.{FlatSpec, Matchers}

class VersionUtilSpec extends FlatSpec with Matchers {

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
