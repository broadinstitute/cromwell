package drs.localizer

import common.util.VersionUtil
import scopt.{OParser, OParserBuilder}


object CommandLineParser {
  val Usage =
    """
    java -jar /path/to/localizer.jar --cloud azure drs://provider/object /local/path/to/file.txt [--vault-name <name>] [--secret-name <name>] [--identity-client-id <id>]
OR
    java -jar /path/to/localizer.jar --cloud google drs://provider/object /local/path/to/file.txt [--requester-pays-project <project>]
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
      opt[String]('c', "cloud").text("Cloud vendor for which to generate an access token").required().
        action((s, c) => c.copy(cloudName = Option(s))),
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
        c.cloudName match {
          case Some("azure") if c.googleRequesterPaysProject.isDefined =>
            failure("'requester-pays-project' is only valid for --cloud google")
          case Some("google") if c.azureVaultName.isDefined =>
            failure("'vault-name' is only valid for --cloud azure")
          case Some("google") if c.azureSecretName.isDefined =>
            failure("'secret-name' is only valid for --cloud azure")
          case Some("google") if c.azureIdentityClientId.isDefined =>
            failure("'identity-client-id' is only valid for --cloud azure")
          case _ => success
        }
      )
    )
  }
}

case class CommandLineArguments(cloudName: Option[String] = None,
                                drsObject: Option[String] = None,
                                containerPath: Option[String] = None,
                                googleRequesterPaysProject: Option[String] = None,
                                azureVaultName: Option[String] = None,
                                azureSecretName: Option[String] = None,
                                azureIdentityClientId: Option[String] = None)
