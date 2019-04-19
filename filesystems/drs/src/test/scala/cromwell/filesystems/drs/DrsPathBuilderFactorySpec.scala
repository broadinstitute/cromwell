package cromwell.filesystems.drs

import com.typesafe.config.ConfigFactory
import cromwell.core.filesystem.CromwellFileSystems
import org.scalatest.{FlatSpec, Matchers}

class DrsPathBuilderFactorySpec extends FlatSpec with Matchers{

  behavior of "DrsPathBuilderFactory"

  it should "create a drs filesystem from a config" in {
    val globalFileSystemConfig = ConfigFactory.parseString(
      """|filesystems {
         |  drs {
         |    class = "cromwell.filesystems.drs.DrsPathBuilderFactory"
         |    global {
         |      class = "cromwell.filesystems.drs.DrsFileSystemConfig"
         |      config {
         |        martha {
         |          url = "http://matha-url"
         |          request.json-template = "{"key": "${holder}"}"
         |        }
         |      }
         |    }
         |  }
         |}
         |""".stripMargin
    )

    val fileSystemConfig = ConfigFactory.parseString(
      """|filesystems {
         |  drs {}
         |}
         |""".stripMargin
    )

    val fileSystems = new CromwellFileSystems(globalFileSystemConfig).factoriesFromConfig(fileSystemConfig).right.get
    fileSystems.keys should contain theSameElementsAs List("drs")
    fileSystems("drs") should be(a[DrsPathBuilderFactory])
  }
}
