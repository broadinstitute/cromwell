package drs.localizer

import common.util.VersionUtil
import drs.localizer.CommandLineParser.AccessTokenStrategy._
import drs.localizer.CommandLineParser.Usage


class CommandLineParser extends scopt.OptionParser[CommandLineArguments](Usage) {
  lazy val localizerVersion: String = VersionUtil.getVersion("cromwell-drs-localizer")

  version("version")

  help("help").text("Cromwell DRS Localizer")

  head("cromwell-drs-localizer", localizerVersion)

  opt[String]('d', "drs-object-id").text("DRS object ID").
    action((s, c) =>
      c.copy(drsObject = Option(s)))
  opt[String]('c',"container-path").text("Container path").
    action((s, c) =>
      c.copy(containerPath = Option(s)))
  opt[String]('m', "manifest-path").text("File path of manifest containing multiple files to localize").
    action((s, c) =>
      c.copy(manifestPath = Option(s)))
  opt[String]('r', "requester-pays-project").text(s"Requester pays project (only valid with '$Google' auth strategy)").optional().
    action((s, c) =>
      c.copy(googleRequesterPaysProject = Option(s)))
  opt[String]('t', "access-token-strategy").text(s"Access token strategy, must be one of '$Azure' or '$Google' (default '$Google')").
    action((s, c) =>
      c.copy(accessTokenStrategy = Option(s.toLowerCase())))
  opt[String]('v', "vault-name").text("Azure vault name").
    action((s, c) =>
      c.copy(azureVaultName = Option(s)))
  opt[String]('s', "secret-name").text("Azure secret name").
    action((s, c) =>
      c.copy(azureSecretName = Option(s)))
  opt[String]('i', "identity-client-id").text("Azure identity client id").
    action((s, c) =>
      c.copy(azureIdentityClientId = Option(s)))
  checkConfig(c =>
    c.accessTokenStrategy match {
      case Some(Azure) if c.googleRequesterPaysProject.nonEmpty =>
        Left(s"Requester pays project is only valid with access token strategy '$Google'")
      case Some(Azure) if List(c.azureVaultName, c.azureSecretName).exists(_.isEmpty) =>
        Left(s"Both vault name and secret name must be specified for access token strategy $Azure")
      case Some(Azure) => Right(())
      case Some(Google) if List(c.azureSecretName, c.azureVaultName, c.azureIdentityClientId).forall(_.isEmpty) => Right(())
      case Some(Google) => Left(s"One or more specified options are only valid with access token strategy '$Azure'")
      case Some(huh) => Left(s"Unrecognized access token strategy '$huh'")
      case None => Left("Programmer error, access token strategy should not be None")
    }
  )
  checkConfig(c =>
    (c.drsObject, c.containerPath, c.manifestPath) match {
      case (Some(_), Some(_), None) => Right(())
      case (None, None, Some(_)) => Right(())
      case _ => Left("Must provide either DRS path and container path, OR manifest file (-m).")
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

    Can be run to localize a single file with DRS id and local container path provided in args:
    java -jar /path/to/localizer.jar [options] -d drs://provider/object -c /local/path/to/file.txt

    Can also be used to localize multiple files in one invocation with manifest file provided in args:
    java -jar /path/to/localizer.jar [options] -m /local/path/to/manifest/file

    Note that the <requester pays project> optional argument is only valid with access token strategy 'Google'.
    """

}

case class CommandLineArguments(accessTokenStrategy: Option[String] = Option(Google),
                                drsObject: Option[String] = None,
                                containerPath: Option[String] = None,
                                googleRequesterPaysProject: Option[String] = None,
                                azureVaultName: Option[String] = None,
                                azureSecretName: Option[String] = None,
                                azureIdentityClientId: Option[String] = None,
                                manifestPath: Option[String] = None)
