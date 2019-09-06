import Testing._
import sbt.Keys._
import sbt._
import sbt.io.Path._

import scala.sys.process._

object ContinuousIntegration {
  val copyCiResources: TaskKey[Unit] = taskKey[Unit](s"Copy CI resources.")
  val renderCiResources: TaskKey[Unit] = taskKey[Unit](s"Render CI resources with Hashicorp Vault.")

  val srcCiResources: SettingKey[File] = settingKey[File]("Source directory for CI resources")
  val targetCiResources: SettingKey[File] = settingKey[File]("Target directory for CI resources")
  val vaultToken: SettingKey[File] = settingKey[File]("File with the vault token")

  val ciSettings: Seq[Setting[_]] = List(
    srcCiResources := sourceDirectory.value / "ci" / "resources",
    targetCiResources := target.value / "ci" / "resources",
    vaultToken := userHome / ".vault-token",
    copyCiResources := {
      IO.copyDirectory(srcCiResources.value, targetCiResources.value)
    },
    renderCiResources := {
      minnieKenny.toTask("").value
      copyCiResources.value
      val log = streams.value.log
      if (!vaultToken.value.exists()) {
        sys.error(
          s"""The vault token file "${vaultToken.value}" does not exist. Be sure to login using the instructions """ +
            """on https://hub.docker.com/r/broadinstitute/dsde-toolbox/ under "Authenticating to vault"."""
        )
      }
      if (vaultToken.value.isDirectory) {
        sys.error(s"""The vault token file "${vaultToken.value}" should not be a directory.""")
      }
      val cmd = List(
        "docker",
        "run",
        "--rm",
        "-v", s"${vaultToken.value}:/root/.vault-token",
        "-v", s"${srcCiResources.value}:${srcCiResources.value}",
        "-v", s"${targetCiResources.value}:${targetCiResources.value}",
        "-e", "ENVIRONMENT=not_used",
        "-e", s"INPUT_PATH=${srcCiResources.value}",
        "-e", s"OUT_PATH=${targetCiResources.value}",
        "broadinstitute/dsde-toolbox:dev", "render-templates.sh"
      )
      val result = cmd ! log
      if (result != 0) {
        sys.error(
          "Vault rendering failed. Please double check for errors above and see the setup instructions on " +
            "https://hub.docker.com/r/broadinstitute/dsde-toolbox/"
        )
      }
    },
  )
}
