package cromwell.engine.io.nio

import akka.actor.ActorRef
import akka.stream.ActorMaterializer
import akka.stream.scaladsl.{Keep, Sink, Source}
import cats.effect.IO
import com.google.cloud.storage.StorageException
import common.mock.MockSugar
import cromwell.core.callcaching.FileHashStrategy
import cromwell.core.io.DefaultIoCommandBuilder._
import cromwell.core.io._
import cromwell.core.path.DefaultPathBuilder
import cromwell.core.{CromwellFatalExceptionMarker, TestKitSuite}
import cromwell.engine.io.IoActor.DefaultCommandContext
import cromwell.engine.io.IoAttempts.EnhancedCromwellIoException
import cromwell.engine.io.IoCommandContext
import cromwell.filesystems.blob.BlobPath
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.http.HttpPathBuilder
import org.mockito.ArgumentMatchers._
import org.mockito.Mockito.when
import org.scalatest.flatspec.AsyncFlatSpecLike
import org.scalatest.matchers.should.Matchers

import java.nio.file.NoSuchFileException
import java.util.UUID
import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.control.NoStackTrace
import scala.util.{Failure, Success, Try}

class NioFlowSpec extends TestKitSuite with AsyncFlatSpecLike with Matchers with MockSugar {

  behavior of "NioFlowSpec"

  private val NoopOnRetry: IoCommandContext[_] => Throwable => Unit = _ => _ => ()
  private val NoopOnBackpressure: Option[Double] => Unit = _ => ()

  private val flow = new NioFlow(parallelism = 1,
                                 onRetryCallback = NoopOnRetry,
                                 onBackpressure = NoopOnBackpressure,
                                 numberOfAttempts = 3,
                                 commandBackpressureStaleness = 5 seconds
  )(system).flow

  implicit val materializer: ActorMaterializer = ActorMaterializer()
  private val replyTo = mock[ActorRef]
  private val readSink = Sink.head[(IoAck[_], IoCommandContext[_])]

  override def afterAll(): Unit = {
    materializer.shutdown()
    super.afterAll()
  }

  it should "write to a Nio Path" in {
    val testPath = DefaultPathBuilder.createTempFile()
    val context = DefaultCommandContext(writeCommand(testPath, "hello", Seq.empty).get, replyTo)
    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)

    stream.run() map { _ =>
      assert(testPath.contentAsString == "hello")
    }
  }

  it should "read from a Nio Path" in {
    val testPath = DefaultPathBuilder.createTempFile()
    testPath.write("hello")

    val context = DefaultCommandContext(contentAsStringCommand(testPath, None, failOnOverflow = true).get, replyTo)
    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)

    stream.run() map {
      case (success: IoSuccess[_], _) => assert(success.result.asInstanceOf[String] == "hello")
      case (ack, _) =>
        fail(s"read returned an unexpected message:\n$ack\n\n")
    }
  }

  it should "get size from a Nio Path" in {
    val testPath = DefaultPathBuilder.createTempFile()
    testPath.write("hello")

    val context = DefaultCommandContext(sizeCommand(testPath).get, replyTo)
    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)

    stream.run() map {
      case (success: IoSuccess[_], _) => assert(success.result.asInstanceOf[Long] == 5)
      case _ => fail("size returned an unexpected message")
    }
  }

  it should "fail with an UnknownHost error when trying to get size for a bogus HTTP path" in {
    val httpPath = new HttpPathBuilder().build("http://ex000mple.c0m/bogus/url/fake.html").get

    val context = DefaultCommandContext(sizeCommand(httpPath).get, replyTo)
    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)
    stream.run() map {
      case (IoFailure(_, EnhancedCromwellIoException(_, receivedException)), _) =>
        receivedException.getMessage should include("UnknownHost")
      case (ack, _) => fail(s"size should have failed with UnknownHost but didn't:\n$ack\n\n")
    }
  }

  it should "fail when trying to get size for a bogus HTTP path" in {
    val httpPath = new HttpPathBuilder().build("http://google.com/bogus/je8934hufe832489uihewuihf").get

    val context = DefaultCommandContext(sizeCommand(httpPath).get, replyTo)
    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)
    stream.run() map {
      case (IoFailure(_, EnhancedCromwellIoException(_, receivedException)), _) =>
        receivedException.getMessage should include("Couldn't fetch size")
      case (ack, _) => fail(s"size should have failed but didn't:\n$ack\n\n")
    }
  }

  it should "get hash from a Nio Path" in {
    val testPath = DefaultPathBuilder.createTempFile()
    testPath.write("hello")

    val context = DefaultCommandContext(hashCommand(testPath, FileHashStrategy.Md5).get, replyTo)
    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)

    stream.run() map {
      case (success: IoSuccess[_], _) =>
        assert(success.result.asInstanceOf[String] == "5d41402abc4b2a76b9719d911017c592")
      case _ => fail("hash returned an unexpected message")
    }
  }

  it should "get hash from a GcsPath" in {
    val exception = new Exception("everything's fine, I am an expected blob failure") with NoStackTrace
    val testPath = mock[GcsPath]
    testPath.objectBlobId returns Failure(exception)

    val context = DefaultCommandContext(hashCommand(testPath, FileHashStrategy.Crc32c).get, replyTo)
    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)

    stream.run() map {
      case (IoFailure(_, EnhancedCromwellIoException(_, receivedException)), _) =>
        receivedException should be(exception)
      case unexpected => fail(s"hash returned an unexpected message: $unexpected")
    }
  }

  it should "get hash from a BlobPath when stored hash exists" in {
    val testPath = mock[BlobPath]
    val hashString = "2d01d5d9c24034d54fe4fba0ede5182d" // echo "hello there" | md5sum
    testPath.md5HexString returns Try(Option(hashString))

    val context = DefaultCommandContext(hashCommand(testPath, FileHashStrategy.Md5).get, replyTo)
    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)

    stream.run() map {
      case (success: IoSuccess[_], _) => assert(success.result.asInstanceOf[String] == hashString)
      case (ack, _) =>
        fail(s"read returned an unexpected message:\n$ack\n\n")
    }
  }

  it should "succeed if a BlobPath is missing a stored hash" in {
    val testPath = mock[BlobPath]
    when(testPath.limitFileContent(any[Option[Int]], any[Boolean])(any[ExecutionContext]))
      .thenReturn("hello there".getBytes)
    when(testPath.md5HexString)
      .thenReturn(Success(None))

    val context =
      DefaultCommandContext(contentAsStringCommand(testPath, Option(100), failOnOverflow = true).get, replyTo)
    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)

    stream.run() map {
      case (success: IoSuccess[_], _) => assert(success.result.asInstanceOf[String] == "hello there")
      case (ack, _) =>
        fail(s"read returned an unexpected message:\n$ack\n\n")
    }
  }

  it should "copy Nio paths" in {
    val testPath = DefaultPathBuilder.createTempFile()
    val testCopyPath = testPath.sibling(UUID.randomUUID().toString)

    val context = DefaultCommandContext(copyCommand(testPath, testCopyPath).get, replyTo)

    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)

    stream.run() map {
      case (_: IoSuccess[_], _) => assert(testCopyPath.exists)
      case _ => fail("copy returned an unexpected message")
    }
  }

  it should "copy Nio paths with" in {
    val testPath = DefaultPathBuilder.createTempFile()
    testPath.write("goodbye")

    val testCopyPath = DefaultPathBuilder.createTempFile()
    testCopyPath.write("hello")

    val context = DefaultCommandContext(copyCommand(testPath, testCopyPath).get, replyTo)

    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)

    stream.run() map {
      case (_: IoSuccess[_], _) =>
        assert(testCopyPath.exists)
        assert(testCopyPath.contentAsString == "goodbye")
      case _ => fail("copy returned an unexpected message")
    }
  }

  it should "delete a Nio path" in {
    val testPath = DefaultPathBuilder.createTempFile()
    val context = DefaultCommandContext(deleteCommand(testPath, swallowIoExceptions = false).get, replyTo)
    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)

    stream.run() map {
      case (_: IoSuccess[_], _) => assert(!testPath.exists)
      case _ => fail("delete returned an unexpected message")
    }
  }

  it should "delete a Nio path with swallowIoExceptions true" in {
    val testPath = DefaultPathBuilder.build("/this/does/not/exist").get

    // noinspection RedundantDefaultArgument
    val context = DefaultCommandContext(deleteCommand(testPath, swallowIoExceptions = true).get, replyTo)

    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)

    stream.run() map {
      case (_: IoSuccess[_], _) => assert(!testPath.exists)
      case _ => fail("delete returned an unexpected message")
    }
  }

  it should "delete a Nio path with swallowIoExceptions false" in {
    val testPath = DefaultPathBuilder.build("/this/does/not/exist").get

    val context = DefaultCommandContext(deleteCommand(testPath, swallowIoExceptions = false).get, replyTo)

    val testSource = Source.single(context)

    val stream = testSource.via(flow).toMat(readSink)(Keep.right)

    stream.run() map {
      case (failure: IoFailure[_], _) =>
        assert(failure.failure.isInstanceOf[CromwellFatalExceptionMarker])
        assert(failure.failure.getMessage == "[Attempted 1 time(s)] - NoSuchFileException: /this/does/not/exist")
        assert(failure.failure.getCause.isInstanceOf[NoSuchFileException])
      case _ => fail(s"delete returned an unexpected message")
    }
  }

  it should "retry on retryable exceptions" in {
    val testPath = DefaultPathBuilder.build("does/not/matter").get

    val context = DefaultCommandContext(contentAsStringCommand(testPath, None, failOnOverflow = false).get, replyTo)

    val testSource = Source.single(context)

    val customFlow = new NioFlow(parallelism = 1,
                                 onRetryCallback = NoopOnRetry,
                                 onBackpressure = NoopOnBackpressure,
                                 numberOfAttempts = 3,
                                 commandBackpressureStaleness = 5 seconds
    )(system) {

      private var tries = 0
      override def handleSingleCommand(ioSingleCommand: IoCommand[_]): IO[IoSuccess[_]] =
        IO {
          tries += 1
          if (tries < 3) throw new StorageException(500, "message")
          else IoSuccess(ioSingleCommand, "content")
        }
    }.flow

    val stream = testSource.via(customFlow).toMat(readSink)(Keep.right)

    stream.run() map {
      case (success: IoSuccess[_], _) => assert(success.result.asInstanceOf[String] == "content")
      case _ => fail("read returned an unexpected message")
    }
  }

}
