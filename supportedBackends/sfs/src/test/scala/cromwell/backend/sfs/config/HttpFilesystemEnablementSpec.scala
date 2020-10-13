package cromwell.backend.sfs.config

import com.typesafe.config.{Config, ConfigFactory}
import common.assertion.CromwellTimeoutSpec
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.sfs.config.HttpFilesystemEnablementSpec._
import cromwell.core.filesystem.CromwellFileSystems
import net.ceedubs.ficus.Ficus._
import org.scalatest.BeforeAndAfterAll
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike

class HttpFilesystemEnablementSpec extends AnyWordSpecLike with CromwellTimeoutSpec with Matchers with BeforeAndAfterAll {
  "The http filesystem on the default Local backend" should {
    "be enabled unless explicitly disabled" in {
      configuredFilesystems(LocalConfig) shouldEqual Set("http")
      configuredFilesystems(LocalConfig.withHttpFilesystem(enabled = true)) shouldEqual Set("http")
      configuredFilesystems(LocalConfig.withHttpFilesystem(enabled = false)) shouldBe empty
    }
  }
}


object HttpFilesystemEnablementSpec {
  def configuredFilesystems(configs: Configurations): Set[String] = {
    val descriptor = new BackendConfigurationDescriptor(configs.local, configs.global) {
      override lazy val cromwellFileSystems = new CromwellFileSystems(configs.global)
    }
    descriptor.configuredPathBuilderFactories.keys.toSet
  }

  val LocalConfig: Configurations = {
    val rawGlobalConfig = ConfigFactory.load()
    // This is only checking the http filesystem, no need to load all those other filesystems with their
    // filesystem classes that live in other subprojects.
    val globalFilesystemsConfig =
    """
      |http {
      |  class = "cromwell.filesystems.http.HttpPathBuilderFactory"
      |}
    """.stripMargin

    val globalFilesystems = ConfigFactory.parseString(globalFilesystemsConfig)
    val updatedGlobalConfig = rawGlobalConfig.withValue("filesystems", globalFilesystems.root())

    val localConfig = updatedGlobalConfig.as[Config]("backend.providers.Local")

    val localFilesystems = ConfigFactory.parseString("http {} ")
    val updatedLocalConfig = localConfig.withValue("filesystems", localFilesystems.root())
    Configurations(local = updatedLocalConfig, global = updatedGlobalConfig)
  }

  case class Configurations(local: Config, global: Config) {
    def withHttpFilesystem(enabled: Boolean): Configurations = {
      val filesystemsConfig =
        s"""
           |http {
           |  enabled: ${enabled.toString}
           |}
        """.stripMargin
      val filesystems = ConfigFactory.parseString(filesystemsConfig)
      val updatedLocalConfig = local.withValue("filesystems", filesystems.root())
      this.copy(global = this.global, local = updatedLocalConfig)
    }
  }
}
