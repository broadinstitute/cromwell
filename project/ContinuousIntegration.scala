import Testing._
import sbt.Keys._
import sbt._
import sbt.io.Path._

import scala.sys.process._

object ContinuousIntegration {
  lazy val ciSettings: Seq[Setting[_]] = List(
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

  def aggregateSettings(rootProject: Project): Seq[Setting[_]] = List(
    // Before compiling, check if the expected projects are aggregated so that they will be compiled-and-tested too.
    compile in Compile := {
      streams.value.log // make sure logger is loaded
      validateAggregatedProjects(rootProject, state.value)
      (compile in Compile).value
    },
  )

  private val copyCiResources: TaskKey[Unit] = taskKey[Unit](s"Copy CI resources.")
  private val renderCiResources: TaskKey[Unit] = taskKey[Unit](s"Render CI resources with Hashicorp Vault.")

  private val srcCiResources: SettingKey[File] = settingKey[File]("Source directory for CI resources")
  private val targetCiResources: SettingKey[File] = settingKey[File]("Target directory for CI resources")
  private val vaultToken: SettingKey[File] = settingKey[File]("File with the vault token")

  /**
    * For "reasons" these projects are excluded from the root aggregation in build.sbt.
    */
  private val unaggregatedProjects = Map(
    "cwlEncoder" -> "not sure what this is",
    "hybridCarboniteMetadataService" -> "not sure why this is excluded",
  )

  /**
    * Get the list of projects defined in build.sbt excluding the passed in root project.
    */
  private def getBuildSbtNames(rootProject: Project, state: State): Set[String] = {
    val extracted = Project.extract(state)
    extracted.structure.units.flatMap({
      case (_, loadedBuildUnit) => loadedBuildUnit.defined.keys
    }).toSet - rootProject.id
  }

  /**
    * Validates that projects are aggregated.
    */
  private def validateAggregatedProjects(rootProject: Project, state: State): Unit = {
    // Get the list of projects explicitly aggregated
    val projectReferences: Seq[ProjectReference] = rootProject.aggregate
    val localProjectReferences = projectReferences collect {
      case localProject: LocalProject => localProject
    }
    val aggregatedNames = localProjectReferences.map(_.project).toSet

    val buildSbtNames = getBuildSbtNames(rootProject, state)
    val missingNames = buildSbtNames.diff(aggregatedNames ++ unaggregatedProjects.keySet).toList.sorted
    if (missingNames.nonEmpty) {
      sys.error(s"There are projects defined in build.sbt that are not aggregated: ${missingNames.mkString(", ")}")
    }

    val falseNames = unaggregatedProjects.filterKeys(aggregatedNames.contains)
    if (falseNames.nonEmpty) {
      val reasons = falseNames.map({case (name, reason) => s"  ${name}: ${reason}"}).mkString("\n")
      sys.error(s"There are projects aggregated in build.sbt that shouldn't be:\n$reasons")
    }
  }
}
