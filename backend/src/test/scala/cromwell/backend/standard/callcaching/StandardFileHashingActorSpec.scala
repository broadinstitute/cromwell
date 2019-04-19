package cromwell.backend.standard.callcaching

import akka.actor.{ActorRef, Props}
import akka.testkit._
import cromwell.backend.standard.callcaching.StandardFileHashingActor.SingleFileHashRequest
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor}
import cromwell.core.TestKitSuite
import cromwell.core.callcaching.HashingFailedMessage
import cromwell.core.io.{IoCommandBuilder, IoHashCommand, PartialIoCommandBuilder}
import cromwell.core.path.{DefaultPathBuilder, Path}
import org.scalatest.{FlatSpecLike, Matchers}
import wom.values.WomSingleFile

import scala.concurrent.TimeoutException
import scala.concurrent.duration._
import scala.util.Try

class StandardFileHashingActorSpec extends TestKitSuite("StandardFileHashingActorSpec") with ImplicitSender
  with FlatSpecLike with Matchers {

  behavior of "StandardFileHashingActor"

  it should "return a failure to the parent when getPath throws an exception" in {
    val parentProbe = TestProbe()
    val params = StandardFileHashingActorSpec.defaultParams()
    val props = Props(new StandardFileHashingActor(params) {
      override def getPath(str: String): Try[Path] = throw new RuntimeException("I am expected during tests")
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
    val parentProbe = TestProbe()
    val params = StandardFileHashingActorSpec.defaultParams()
    val props = Props(new StandardFileHashingActor(params) {
      override val ioCommandBuilder = IoCommandBuilder(
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
    val parentProbe = TestProbe()
    val ioActorProbe = TestProbe()
    val params = StandardFileHashingActorSpec.ioActorParams(ioActorProbe.ref)
    val props = Props(new StandardFileHashingActor(params) {
      override lazy val defaultIoTimeout = 1.second.dilated

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

}

object StandardFileHashingActorSpec {
  private def testing: Nothing = throw new UnsupportedOperationException("should not be run during tests")

  def defaultParams(): StandardFileHashingActorParams = defaultParams(testing, testing, testing, testing, testing)

  def ioActorParams(ioActor: ActorRef): StandardFileHashingActorParams = {
    defaultParams(withJobDescriptor = testing,
      withConfigurationDescriptor = testing,
      withIoActor = ioActor,
      withServiceRegistryActor = testing,
      withBackendInitializationDataOption = testing)
  }

  def defaultParams(withJobDescriptor: => BackendJobDescriptor,
                    withConfigurationDescriptor: => BackendConfigurationDescriptor,
                    withIoActor: => ActorRef,
                    withServiceRegistryActor: => ActorRef,
                    withBackendInitializationDataOption: => Option[BackendInitializationData]
                   ): StandardFileHashingActorParams = new StandardFileHashingActorParams {

    override def jobDescriptor = withJobDescriptor

    override def configurationDescriptor = withConfigurationDescriptor

    override def ioActor = withIoActor

    override def serviceRegistryActor = withServiceRegistryActor

    override def backendInitializationDataOption = withBackendInitializationDataOption

    override def fileHashCachingActor: Option[ActorRef] = None
  }

}

