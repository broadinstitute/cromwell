package drs.localizer

import common.util.VersionUtil
import drs.localizer.CommandLineParser.AccessTokenStrategy._
import scopt.{OParser, OParserBuilder}


object CommandLineParser {
  object AccessTokenStrategy {
    val Azure = "azure"
    val Google = "google"
  }

  val Usage =
    s"""
    java -jar /path/to/localizer.jar --token-strategy $Azure drs://provider/object /local/path/to/file.txt [--vault-name <name>] [--secret-name <name>] [--identity-client-id <id>]
OR
    java -jar /path/to/localizer.jar --token-strategy $Google drs://provider/object /local/path/to/file.txt [--requester-pays-project <project>]
    """

  lazy val localizerVersion: String = VersionUtil.getVersion("cromwell-drs-localizer")

  def buildParser(): OParser[Unit, CommandLineArguments] = {
    import scopt.OParser
    val builder: OParserBuilder[CommandLineArguments] = OParser.builder[CommandLineArguments]

    import builder._
    OParser.sequence(
      version("version"),
      help("help").text("Cromwell DRS Localizer"),
      head("cromwell-drs-localizer", localizerVersion),
      arg[String]("drs-object-id").text("DRS object ID").required().
        action((s, c) =>
          c.copy(drsObject = Option(s))),
      arg[String]("container-path").text("Container path").required().
        action((s, c) =>
          c.copy(containerPath = Option(s))),
      opt[String]('t', "token-strategy").text("Strategy to use when generating an access token to call Martha").required().
        action((s, c) => c.copy(tokenStrategy = Option(s.toLowerCase()))),
      opt[String]('v', "vault-name").text("Azure vault name").
        action((s, c) =>
          c.copy(azureVaultName = Option(s))),
      opt[String]('s', "secret-name").text("Azure secret name").
        action((s, c) =>
          c.copy(azureSecretName = Option(s))),
      opt[String]('i', "identity-client-id").text("Azure identity client id").
        action((s, c) =>
          c.copy(azureIdentityClientId = Option(s))),
      opt[String]('r', "requester-pays-project").text("Google requester pays project name").
        action((s, c) =>
          c.copy(googleRequesterPaysProject = Option(s))),
      checkConfig(c =>
        c.tokenStrategy match {
          case Some(Azure) if c.googleRequesterPaysProject.isDefined =>
            failure(s"'requester-pays-project' is only valid for --token-strategy $Google")
          case Some(Google) if c.azureVaultName.isDefined || c.azureSecretName.isDefined || c.azureIdentityClientId.isDefined =>
            failure(s"token-strategy $Google specified, but 'vault-name', 'secret-name', and 'identity-client-id' are valid only with token-strategy $Azure")
          case Some(Google) | Some(Azure) => success
          case Some(other) => failure(s"Unrecognized token-strategy '$other', only '$Azure' and '$Google' are supported.")
          case _ => failure("")
        }
      )
    )
  }
}

case class CommandLineArguments(tokenStrategy: Option[String] = None,
                                drsObject: Option[String] = None,
                                containerPath: Option[String] = None,
                                googleRequesterPaysProject: Option[String] = None,
                                azureVaultName: Option[String] = None,
                                azureSecretName: Option[String] = None,
                                azureIdentityClientId: Option[String] = None)
