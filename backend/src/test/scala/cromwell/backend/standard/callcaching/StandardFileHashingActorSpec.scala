package cromwell.backend.standard.callcaching

import akka.actor.{ActorRef, Props}
import akka.testkit._
import com.typesafe.config.ConfigFactory
import cromwell.backend.standard.callcaching.StandardFileHashingActor.SingleFileHashRequest
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor}
import cromwell.core.TestKitSuite
import cromwell.core.callcaching.{FileHashStrategy, HashingFailedMessage, HashType, SuccessfulHashResultMessage}
import cromwell.core.io.{IoCommand, IoCommandBuilder, IoHashCommand, IoSuccess, PartialIoCommandBuilder}
import cromwell.core.path.{DefaultPathBuilder, Path}
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import common.mock.MockSugar
import wom.values.WomSingleFile

import scala.concurrent.TimeoutException
import scala.concurrent.duration._
import scala.util.control.NoStackTrace
import scala.util.{Failure, Success, Try}

class StandardFileHashingActorSpec
    extends TestKitSuite
    with ImplicitSender
    with AnyFlatSpecLike
    with Matchers
    with MockSugar {

  behavior of "StandardFileHashingActor"

  it should "return a failure to the parent when getPath returns an exception" in {
    val parentProbe = TestProbe("parentProbe")
    val params = StandardFileHashingActorSpec.defaultParams()
    val props = Props(new StandardFileHashingActor(params) {
      override def getPath(str: String): Try[Path] =
        Failure(new RuntimeException("I am expected during tests") with NoStackTrace)
    })
    val standardFileHashingActorRef = TestActorRef(props, parentProbe.ref)
    val request = SingleFileHashRequest(null, null, WomSingleFile("/expected/failure/path"), None)
    standardFileHashingActorRef ! request

    parentProbe.expectMsgPF(1.seconds) {
      case failed: HashingFailedMessage if failed.file == "/expected/failure/path" =>
        failed.reason should be(a[RuntimeException])
        failed.reason.getMessage should be("I am expected during tests")
      case unexpected =>
        fail(s"received unexpected message $unexpected")
    }
  }

  it should "return a failure to the parent when hashCommand throws an exception" in {
    val parentProbe = TestProbe("parentProbe")
    val params = StandardFileHashingActorSpec.defaultParams()
    val props = Props(new StandardFileHashingActor(params) {
      override val ioCommandBuilder: IoCommandBuilder = IoCommandBuilder(
        new PartialIoCommandBuilder {
          override def hashCommand = throw new RuntimeException("I am expected during tests")
        }
      )
      override def getPath(str: String): Try[Path] = Try(DefaultPathBuilder.get(str))
    })
    val standardFileHashingActorRef = TestActorRef(props, parentProbe.ref)
    val request = SingleFileHashRequest(null, null, WomSingleFile("/expected/failure/path"), None)
    standardFileHashingActorRef ! request

    parentProbe.expectMsgPF(10.seconds.dilated) {
      case failed: HashingFailedMessage if failed.file == "/expected/failure/path" =>
        failed.reason should be(a[RuntimeException])
        failed.reason.getMessage should be("I am expected during tests")
      case unexpected =>
        fail(s"received unexpected message $unexpected")
    }
  }

  it should "send a timeout to the ioActor the command doesn't hash" in {
    val parentProbe = TestProbe("parentProbe")
    val ioActorProbe = TestProbe("ioActorProbe")
    val params = StandardFileHashingActorSpec.ioActorParams(ioActorProbe.ref)
    val props = Props(new StandardFileHashingActor(params) {
      override lazy val defaultIoTimeout: FiniteDuration = 1.second.dilated

      override def getPath(str: String): Try[Path] = Try(DefaultPathBuilder.get(str))
    })
    val standardFileHashingActorRef = TestActorRef(props, parentProbe.ref)
    val request = SingleFileHashRequest(null, null, WomSingleFile("/expected/failure/path"), None)

    standardFileHashingActorRef ! request

    ioActorProbe.expectMsgPF(10.seconds.dilated) {
      case (request: FileHashContext, _: IoHashCommand) if request.file == "/expected/failure/path" =>
      case unexpected => fail(s"received unexpected message $unexpected")
    }

    parentProbe.expectMsgPF(10.seconds.dilated) {
      case failed: HashingFailedMessage if failed.file == "/expected/failure/path" =>
        failed.reason should be(a[TimeoutException])
        failed.reason.getMessage should be("Hashing request timed out for: /expected/failure/path")
      case unexpected => fail(s"received unexpected message $unexpected")
    }
  }

  it should "handle non-string hash responses" in {
    val parentProbe = TestProbe("testParentHashString")
    val params = StandardFileHashingActorSpec.ioActorParams(ActorRef.noSender)
    val props = Props(new StandardFileHashingActor(params) {
      override lazy val defaultIoTimeout: FiniteDuration = 1.second.dilated

      override def getPath(str: String): Try[Path] = Try(DefaultPathBuilder.get(str))
    })
    val standardFileHashingActorRef = parentProbe.childActorOf(props, "testStandardFileHashingActorHashString")

    val fileHashContext = mock[FileHashContext]
    fileHashContext.file returns "/expected/failure/path"
    val command = mock[IoCommand[Int]]
    val message: (FileHashContext, IoSuccess[Int]) = (fileHashContext, IoSuccess(command, 1357))

    standardFileHashingActorRef ! message

    parentProbe.expectMsgPF(10.seconds.dilated) {
      case failed: HashingFailedMessage =>
        failed.reason should be(a[Exception])
        failed.reason.getMessage should
          be("Hash function supposedly succeeded but responded with '1357' instead of a string hash")
      case unexpected => fail(s"received unexpected message $unexpected")
    }
  }

  it should "handle string hash responses" in {
    val parentProbe = TestProbe("testParentHashString")
    val params = StandardFileHashingActorSpec.ioActorParams(ActorRef.noSender)
    val props = Props(new StandardFileHashingActor(params) {
      override lazy val defaultIoTimeout: FiniteDuration = 1.second.dilated

      override def getPath(str: String): Try[Path] = Try(DefaultPathBuilder.get(str))
    })
    val standardFileHashingActorRef = parentProbe.childActorOf(props, "testStandardFileHashingActorHashString")

    val fileHashContext = mock[FileHashContext]
    fileHashContext.file returns "/expected/failure/path"
    val command = mock[IoCommand[String]]
    val message: (FileHashContext, IoSuccess[String]) = (fileHashContext, IoSuccess(command, "a_nice_hash"))

    standardFileHashingActorRef ! message

    parentProbe.expectMsgPF(10.seconds.dilated) {
      case succeeded: SuccessfulHashResultMessage =>
        succeeded.hashes.map(_.hashValue.value).headOption shouldBe Some("a_nice_hash")
      case unexpected => fail(s"received unexpected message $unexpected")
    }
  }

  it should "use the right hashing strategies" in {
    val parentProbe = TestProbe("testParentHashStrategies")
    val ioActorProbe = TestProbe("ioActorProbe")
    val backendConfig = ConfigFactory.parseString(
      """filesystems.gcs.caching.hashing-strategy = ["md5", "identity"]
        |filesystems.s3.caching.hashing-strategy = "etag"
        |filesystems.http.some-other-config = "foobar"
        |filesystems.ftp.caching.hashing-strategy = []
        |filesystems.nfs.caching.hashing-strategy = "bogohash"""".stripMargin
    )
    val config = BackendConfigurationDescriptor(backendConfig, ConfigFactory.empty)

    val props =
      Props(new StandardFileHashingActor(StandardFileHashingActorSpec.ioActorParams(ioActorProbe.ref, config)) {
        override val defaultHashingStrategies: Map[String, FileHashStrategy] = Map(
          "gcs" -> FileHashStrategy.Crc32c,
          "drs" -> FileHashStrategy.Drs
        )
        override val fallbackHashingStrategy: FileHashStrategy = FileHashStrategy(List(HashType.Sha256))

        override def getPath(str: String): Try[Path] = {
          val p = mock[Path]
          p.filesystemTypeKey returns str
          Success(p)
        }
      })
    val standardFileHashingActorRef = parentProbe.childActorOf(props, "testStandardFileHashingActorHashStrategy")

    def checkHashStrategy(filesystemKey: String, expectedStrategy: FileHashStrategy): Unit = {
      val request = SingleFileHashRequest(null, null, WomSingleFile(filesystemKey), None)
      standardFileHashingActorRef ! request
      ioActorProbe.expectMsgPF(10.seconds.dilated) {
        case (_: FileHashContext, cmd: IoHashCommand) if cmd.hashStrategy == expectedStrategy =>
        case unexpected => fail(s"received unexpected ${filesystemKey} message $unexpected")
      }
    }

    // Test an actor-defined default overriden by config
    checkHashStrategy("gcs", FileHashStrategy(List(HashType.Md5, HashType.Identity)))

    // Test a strategy only defined in config
    checkHashStrategy("s3", FileHashStrategy.ETag)

    // Test a strategy defined as an empty list in config, should use fallback
    checkHashStrategy("ftp", FileHashStrategy(List(HashType.Sha256)))

    // Test a strategy with a default defined in the actor
    checkHashStrategy("drs", FileHashStrategy.Drs)

    // Test a filesystem that has config, but not this config
    checkHashStrategy("http", FileHashStrategy(List(HashType.Sha256)))

    // Test a strategy not defined in config or actor defaults
    checkHashStrategy("blob", FileHashStrategy(List(HashType.Sha256)))

    // Test a strategy with an invalid value in config
    checkHashStrategy("nfs", FileHashStrategy(List(HashType.Sha256)))
  }

}

object StandardFileHashingActorSpec {
  private def testing: Nothing = throw new UnsupportedOperationException("should not be run during tests")
  private val emptyBackendConfig = BackendConfigurationDescriptor(ConfigFactory.empty, ConfigFactory.empty)

  def defaultParams(config: BackendConfigurationDescriptor = emptyBackendConfig): StandardFileHashingActorParams =
    defaultParams(testing, emptyBackendConfig, testing, testing, testing)

  def ioActorParams(ioActor: ActorRef,
                    config: BackendConfigurationDescriptor = emptyBackendConfig
  ): StandardFileHashingActorParams =
    defaultParams(
      withJobDescriptor = testing,
      withConfigurationDescriptor = config,
      withIoActor = ioActor,
      withServiceRegistryActor = testing,
      withBackendInitializationDataOption = testing
    )

  def defaultParams(withJobDescriptor: => BackendJobDescriptor,
                    withConfigurationDescriptor: => BackendConfigurationDescriptor,
                    withIoActor: => ActorRef,
                    withServiceRegistryActor: => ActorRef,
                    withBackendInitializationDataOption: => Option[BackendInitializationData]
  ): StandardFileHashingActorParams = new StandardFileHashingActorParams {

    override def jobDescriptor: BackendJobDescriptor = withJobDescriptor

    override def configurationDescriptor: BackendConfigurationDescriptor = withConfigurationDescriptor

    override def ioActor: ActorRef = withIoActor

    override def serviceRegistryActor: ActorRef = withServiceRegistryActor

    override def backendInitializationDataOption: Option[BackendInitializationData] =
      withBackendInitializationDataOption

    override def fileHashCachingActor: Option[ActorRef] = None
  }

  class StrategyTestFileHashingActor(standardParams: StandardFileHashingActorParams)
      extends StandardFileHashingActor(standardParams) {
    override val defaultHashingStrategies: Map[String, FileHashStrategy] = Map(
      "gcs" -> FileHashStrategy.Crc32c,
      "drs" -> FileHashStrategy.Drs
    )
    override val fallbackHashingStrategy: FileHashStrategy = FileHashStrategy(List(HashType.Sha256))
  }

}
