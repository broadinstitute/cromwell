package cromwell.core.filesystem

import akka.actor.ActorSystem
import cats.data.NonEmptyList
import com.typesafe.config.{Config, ConfigFactory}
import common.exception.AggregatedMessageException
import cromwell.core.WorkflowOptions
import cromwell.core.path.MockPathBuilderFactory
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.ExecutionContext

class CromwellFileSystemsSpec extends FlatSpec with Matchers {
  behavior of "CromwellFileSystems"

  val globalConfig = ConfigFactory.parseString(
    """
      |filesystems {
      |  fs1.class = "cromwell.core.path.MockPathBuilderFactory"
      |  fs2.class = "cromwell.core.path.MockPathBuilderFactory"
      |  fs3.class = "cromwell.core.filesystem.MockNotPathBuilderFactory"
      |}
    """.stripMargin)

  val cromwellFileSystems = new CromwellFileSystems(globalConfig)
  
  it should "build factory builders and factories for valid configuration" in {
    cromwellFileSystems.factoryBuilders.keySet shouldBe Set("fs1", "fs2", "fs3")
    
    val factoriesConfig = ConfigFactory.parseString(
      """
        |filesystems {
        | fs1.somekey = "somevalue"
        | fs2.someotherkey = "someothervalue"
        |}
      """.stripMargin)

    val pathFactories = cromwellFileSystems.factoriesFromConfig(factoriesConfig)
    pathFactories.isRight shouldBe true
    val fs1 = pathFactories.right.get("fs1")
    val fs2 = pathFactories.right.get("fs2")
    fs1 shouldBe a[MockPathBuilderFactory]
    fs2 shouldBe a[MockPathBuilderFactory]
    fs1.asInstanceOf[MockPathBuilderFactory].instanceConfig.getString("somekey") shouldBe "somevalue"
    fs2.asInstanceOf[MockPathBuilderFactory].instanceConfig.getString("someotherkey") shouldBe "someothervalue"
  }
  
  it should "build singleton instance if specified" in {
    val rootConf = ConfigFactory.parseString(
      """
        |filesystems {
        |  fs1 {
        |    class = "cromwell.core.filesystem.MockPathBuilderFactoryCustomSingletonConfig"
        |    global {
        |      class = "cromwell.core.filesystem.MockSingletonConfig"
        |    }
        |  }
        |}
      """.stripMargin)

    val cromwellFileSystems = new CromwellFileSystems(rootConf)
    val constructorAndSingleton = cromwellFileSystems.factoryBuilders("fs1")
    constructorAndSingleton._2.isDefined shouldBe true
    constructorAndSingleton._2.get.isInstanceOf[MockSingletonConfig] shouldBe true

    val factory1 = cromwellFileSystems.buildFactory("fs1", ConfigFactory.empty)
    val factory2 = cromwellFileSystems.buildFactory("fs1", ConfigFactory.empty)

    // The singleton configs should be the same for different factories
    assert(factory1.right.get.asInstanceOf[MockPathBuilderFactoryCustomSingletonConfig].singletonConfig ==
      factory2.right.get.asInstanceOf[MockPathBuilderFactoryCustomSingletonConfig].singletonConfig)
  }

  List(
    ("if the filesystem does not exist", "filesystems.fs4.key = value", NonEmptyList.one("Cannot find a filesystem with name fs4 in the configuration. Available filesystems: fs1, fs2, fs3")),
    ("if the config is invalid", "filesystems.fs1 = true", NonEmptyList.one("Invalid filesystem backend configuration for fs1")),
    ("the class is not a PathBuilderFactory", "filesystems.fs3.key = value", NonEmptyList.one("The filesystem class for fs3 is not an instance of PathBuilderFactory"))
  ) foreach {
    case (description, config, expected) =>
      it should s"fail to build factories $description" in {
        val result = cromwellFileSystems.factoriesFromConfig(ConfigFactory.parseString(config))
          result.isLeft shouldBe true
        result.left.get shouldBe expected
      }
  }
  
  val classNotFoundException = AggregatedMessageException(
    "Failed to initialize Cromwell filesystems",
    List("Class do.not.exists for filesystem fs1 cannot be found in the class path.")
  )

  val wrongSignatureException = AggregatedMessageException(
    "Failed to initialize Cromwell filesystems",
    List("Class cromwell.core.filesystem.MockPathBuilderFactoryWrongSignature for filesystem fs1 does not have the required constructor signature: (com.typesafe.config.Config, com.typesafe.config.Config)")
  )

  val invalidConfigException = AggregatedMessageException(
    "Failed to initialize Cromwell filesystems",
    List("Invalid filesystem configuration for gcs")
  )

  val missingClassFieldException = AggregatedMessageException(
    "Failed to initialize Cromwell filesystems",
    List("Filesystem configuration fs1 doesn't have a class field")
  )
  
  List(
    ("is invalid", "filesystems.gcs = true", invalidConfigException),
    ("is missing class fields", "filesystems.fs1.notclass = hello", missingClassFieldException),
    ("can't find class", "filesystems.fs1.class = do.not.exists", classNotFoundException),
    ("has invalid class signature", "filesystems.fs1.class = cromwell.core.filesystem.MockPathBuilderFactoryWrongSignature", wrongSignatureException)
  ) foreach {
    case (description, config, expected) =>
      it should s"fail if global filesystems config $description" in {
        val ex = the[Exception] thrownBy { new CromwellFileSystems(ConfigFactory.parseString(config)) }
        ex shouldBe expected
      }
  }
}

class MockPathBuilderFactoryWrongSignature()
class MockNotPathBuilderFactory(globalConfig: Config, val instanceConfig: Config)

class MockSingletonConfig(config: Config)
class MockPathBuilderFactoryCustomSingletonConfig(globalConfig: Config, val instanceConfig: Config, val singletonConfig: MockSingletonConfig) extends cromwell.core.path.PathBuilderFactory {
  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext) = throw new UnsupportedOperationException
}
