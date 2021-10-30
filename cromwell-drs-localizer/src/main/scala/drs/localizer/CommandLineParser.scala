package drs.localizer

import common.util.VersionUtil


object CommandLineParser {
  val Usage =
    """
    java -jar /path/to/localizer.jar --cloud azure drs://provider/object /local/path/to/file.txt [--vault-name <name>] [--secret-name <name>] [--identity-client-id <id>]
OR
    java -jar /path/to/localizer.jar --cloud google drs://provider/object /local/path/to/file.txt [--requester-pays-project <project>]
    """

  lazy val localizerVersion: String = VersionUtil.getVersion("cromwell-drs-localizer")
}

case class CommandLineArguments(cloudName: Option[String] = None,
                                drsObject: Option[String] = None,
                                containerPath: Option[String] = None,
                                googleRequesterPaysProject: Option[String] = None,
                                azureVaultName: Option[String] = None,
                                azureSecretName: Option[String] = None,
                                azureIdentityClientId: Option[String] = None)
