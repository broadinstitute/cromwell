package drs.localizer

import common.util.VersionUtil
import drs.localizer.CommandLineParser.localizerVersion
import drs.localizer.CommandLineParser.Usage


object CommandLineParser {
  val Usage =
    """
    java -jar /path/to/localizer.jar azure drs://provider/object /local/path/to/file.txt [--vault-name <name>] [--secret-name <name>] [--identity-client-id <id>]
OR
    java -jar /path/to/localizer.jar google drs://provider/object /local/path/to/file.txt [--requester-pays-project <project>]
    """

  lazy val localizerVersion: String = VersionUtil.getVersion("cromwell-drs-localizer")
}

class CommandLineParser extends scopt.OptionParser[CommandLineArguments](Usage) {
  def commonArguments = List(
    arg[String]("drs-id").text("DRS object ID").required().
      action((s, c) =>
        c.copy(drsObject = Option(s))),
    arg[String]("container-path").text("Container path").required().
      action((s, c) =>
        c.copy(containerPath = Option(s))),
  )

  head("cromwell-drs-localizer", localizerVersion)

  help("help").text("Cromwell DRS Localizer")

  version("version")

  cmd("azure").action((_, c) => c.copy(cloudName = Option("azure"))).
    text("Generates an access token for an Azure Batch VM environment").
    children(
      commonArguments ++ List(
        opt[String]('v', "vault-name").text("Azure vault name").
          action((s, c) =>
            c.copy(azureVaultName = Option(s))),
        opt[String]('s', "secret-name").text("Azure secret name").
          action((s, c) =>
            c.copy(azureSecretName = Option(s))),
        opt[String]('i', "identity-client-id").text("Azure identity client id").
          action((s, c) =>
            c.copy(azureIdentityClientId = Option(s)))
      ): _*
    )

  cmd("google").action((_, c) => c.copy(cloudName = Option("google"))).
    text("Generates an access token for a Google Genomics / Lifesciences VM environment").
    children(
      commonArguments ++ List(
        opt[String]('r', "requester-pays-project").text("Google requester pays project name").
          action((s, c) =>
            c.copy(googleRequesterPaysProject = Option(s)))
      ): _*
    )
}

case class CommandLineArguments(cloudName: Option[String] = None,
                                drsObject: Option[String] = None,
                                containerPath: Option[String] = None,
                                googleRequesterPaysProject: Option[String] = None,
                                azureVaultName: Option[String] = None,
                                azureSecretName: Option[String] = None,
                                azureIdentityClientId: Option[String] = None)
