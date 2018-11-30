//package cromwell.filesystems.demo.dos
//
//import com.typesafe.config.ConfigFactory
//import cromwell.core.filesystem.CromwellFileSystems
//import org.scalatest.{FlatSpec, Matchers}
//
//class DemoDosPathBuilderFactorySpec extends FlatSpec with Matchers {
//
//  behavior of "DemoDosPathBuilderFactory"
//
//  it should "create a demo-dos filesystem from a config" in {
//    val globalConfig = ConfigFactory.parseString(
//      """|filesystems {
//         |  demo-dos {
//         |    class = "cromwell.filesystems.demo.dos.DemoDosPathBuilderFactory"
//         |  }
//         |}
//         |""".stripMargin
//    )
//    val fileSystemConfig = ConfigFactory.parseString(
//      """|filesystems {
//         |  demo-dos {}
//         |}
//         |""".stripMargin
//    )
//    val fileSystems = new CromwellFileSystems(globalConfig).factoriesFromConfig(fileSystemConfig).right.get
//    fileSystems.keys should contain theSameElementsAs List("demo-dos")
//    fileSystems("demo-dos") should be(a[DemoDosPathBuilderFactory])
//  }
//
//}
