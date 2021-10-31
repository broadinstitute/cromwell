package drs.localizer

import common.util.VersionUtil
import drs.localizer.CommandLineParser.AccessTokenStrategy._
import drs.localizer.CommandLineParser.Usage


class CommandLineParser extends scopt.OptionParser[CommandLineArguments](Usage) {
  lazy val localizerVersion: String = VersionUtil.getVersion("cromwell-drs-localizer")

  val commonArguments = List(
    arg[String]("drs-object-id").text("DRS object ID").required().
      action((s, c) =>
        c.copy(drsObject = Option(s))),
    arg[String]("container-path").text("Container path").required().
      action((s, c) =>
        c.copy(containerPath = Option(s))),
  )

  version("version")

  help("help").text("Cromwell DRS Localizer")

  head("cromwell-drs-localizer", localizerVersion)

  cmd("azure").
    text("Localize DRS file using Azure UAMI / B2C access token strategy").
    children(commonArguments ++ List(
      opt[String]('v', "vault-name").text("Azure vault name").
        action((s, c) =>
          c.copy(azureVaultName = Option(s))),
      opt[String]('s', "secret-name").text("Azure secret name").
        action((s, c) =>
          c.copy(azureSecretName = Option(s))),
      opt[String]('i', "identity-client-id").text("Azure identity client id").
        action((s, c) =>
          c.copy(azureIdentityClientId = Option(s))),
    ): _*)

  cmd("google").
    text("Localize DRS file using Google Application Default Credentials access token strategy").
    children(commonArguments ++ List(
      opt[String]('r', "requester-pays-project").text("Google requester pays project name").
        action((s, c) =>
          c.copy(googleRequesterPaysProject = Option(s))),
    ): _*)
}

object CommandLineParser {
  /**
   * These access token strategies are named simplistically as there is currently only one access token strategy being
   * used for each of these cloud vendors. But it is certainly possible that multiple strategies could come into use
   * for a particular vendor, in which case the names may need to become more specific for disambiguation.
   */
  object AccessTokenStrategy {
    val Azure = "azure"
    val Google = "google"
  }

  val Usage =
    s"""
Usage:
    java -jar /path/to/localizer.jar $Azure drs://provider/object /local/path/to/file.txt [--vault-name <name>] [--secret-name <name>] [--identity-client-id <id>]
OR
    java -jar /path/to/localizer.jar $Google drs://provider/object /local/path/to/file.txt [--requester-pays-project <project>]
    """

}

case class CommandLineArguments(accessTokenStrategy: Option[String] = None,
                                drsObject: Option[String] = None,
                                containerPath: Option[String] = None,
                                googleRequesterPaysProject: Option[String] = None,
                                azureVaultName: Option[String] = None,
                                azureSecretName: Option[String] = None,
                                azureIdentityClientId: Option[String] = None)
