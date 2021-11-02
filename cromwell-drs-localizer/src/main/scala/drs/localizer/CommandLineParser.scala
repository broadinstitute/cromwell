package drs.localizer

import common.util.VersionUtil
import drs.localizer.CommandLineParser.AccessTokenStrategy._
import drs.localizer.CommandLineParser.Usage


class CommandLineParser extends scopt.OptionParser[CommandLineArguments](Usage) {
  lazy val localizerVersion: String = VersionUtil.getVersion("cromwell-drs-localizer")

  version("version")

  help("help").text("Cromwell DRS Localizer")

  head("cromwell-drs-localizer", localizerVersion)

  arg[String]("drs-object-id").text("DRS object ID").required().
    action((s, c) =>
      c.copy(drsObject = Option(s)))
  arg[String]("container-path").text("Container path").required().
    action((s, c) =>
      c.copy(containerPath = Option(s)))
  opt[String]('t', "access-token-strategy").text(s"Access token strategy, must be one of '$Azure' or '$Google''").
    action((s, c) =>
      c.copy(accessTokenStrategy = Option(s)))
  opt[String]('v', "vault-name").text("Azure vault name").
    action((s, c) =>
      c.copy(azureVaultName = Option(s)))
  opt[String]('s', "secret-name").text("Azure secret name").
    action((s, c) =>
      c.copy(azureSecretName = Option(s)))
  opt[String]('i', "identity-client-id").text("Azure identity client id").
    action((s, c) =>
      c.copy(azureIdentityClientId = Option(s)))
  opt[String]('r', "requester-pays-project").text("Google requester pays project name").
    action((s, c) =>
      c.copy(googleRequesterPaysProject = Option(s)))
  checkConfig(c =>
    c.accessTokenStrategy match {
      case Some(Azure) if c.googleRequesterPaysProject.isEmpty => Right(())
      case Some(Google) if List(c.azureSecretName, c.azureVaultName, c.azureIdentityClientId).forall(_.isEmpty) => Right(())
      case Some(Azure) => Left(s"Requester pays project is only valid with access token strategy $Google")
      case Some(Google) => Left(s"One or more specified options are only valid with access token strategy $Azure")
      case Some(huh) => Left(s"Unrecognized access token strategy $huh")
      case None => Left("Unspecified access token strategy")
    }
  )
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
    java -jar /path/to/localizer.jar [options] drs://provider/object /local/path/to/file.txt
    """

}

case class CommandLineArguments(accessTokenStrategy: Option[String] = Option(Google),
                                drsObject: Option[String] = None,
                                containerPath: Option[String] = None,
                                googleRequesterPaysProject: Option[String] = None,
                                azureVaultName: Option[String] = None,
                                azureSecretName: Option[String] = None,
                                azureIdentityClientId: Option[String] = None)
