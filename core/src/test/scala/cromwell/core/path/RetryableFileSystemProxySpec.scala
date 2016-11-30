package cromwell.core.path

import java.io.FileNotFoundException
import java.nio.channels.SeekableByteChannel
import java.nio.file.DirectoryStream.Filter
import java.nio.file.attribute.{BasicFileAttributes, FileAttributeView}
import java.nio.file.spi.FileSystemProvider
import java.nio.file.{DirectoryStream, OpenOption, Path, StandardOpenOption}
import java.util.concurrent.TimeoutException

import cromwell.core.path.proxy.RetryableFileSystemProviderProxy
import cromwell.core.retry.Backoff
import cromwell.core.{CromwellFatalException, TestKitSuite}
import org.mockito.Matchers._
import org.mockito.Mockito._
import org.mockito.invocation.InvocationOnMock
import org.mockito.stubbing.Answer
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps

class RetryableFileSystemProxySpec extends TestKitSuite with FlatSpecLike with Matchers {

  behavior of "RetryableFileSystemProxySpec"

  case class ThrowParams(exception: Exception, nbTimes: Int)

  abstract class FileSystemAnswer[T](delay: Option[Duration] = None,
                                     throws: Option[ThrowParams] = None) extends Answer[T] {

    var nbThrows = 0

    def delayAndOrThrow() = {
      delay foreach { d => Thread.sleep(d.toMillis) }
      throws foreach { e =>
        if (nbThrows < e.nbTimes) {
          nbThrows = nbThrows + 1
          throw e.exception
        }
      }
    }
  }

  def mockFileSystem(delay: Option[Duration] = None,
                     throws: Option[ThrowParams] = None): FileSystemProvider =  {

    val provider = mock(classOf[FileSystemProvider])

    def answerUnit: Answer[Unit] = new FileSystemAnswer[Unit](delay, throws) {
      override def answer(invocation: InvocationOnMock): Unit = delayAndOrThrow()
    }

    def answerBoolean: Answer[Boolean] = new FileSystemAnswer[Boolean](delay, throws) {
      override def answer(invocation: InvocationOnMock): Boolean = {
        delayAndOrThrow()
        true
      }
    }

    def answerSeekableByteChannel: Answer[SeekableByteChannel] = new FileSystemAnswer[SeekableByteChannel](delay, throws) {
      override def answer(invocation: InvocationOnMock): SeekableByteChannel = {
        delayAndOrThrow()
        mock(classOf[SeekableByteChannel])
      }
    }

    def answerDirectoryStream: Answer[DirectoryStream[Path]] = new FileSystemAnswer[DirectoryStream[Path]](delay, throws) {
      override def answer(invocation: InvocationOnMock): DirectoryStream[Path] = {
        delayAndOrThrow()
        mock(classOf[DirectoryStream[Path]])
      }
    }

    def answerBasicFileAttributes: Answer[BasicFileAttributes] = new FileSystemAnswer[BasicFileAttributes](delay, throws) {
      override def answer(invocation: InvocationOnMock): BasicFileAttributes = {
        delayAndOrThrow()
        mock(classOf[BasicFileAttributes])
      }
    }

    def answerMap: Answer[java.util.Map[String, AnyRef]] = new FileSystemAnswer[java.util.Map[String, AnyRef]](delay, throws) {
      override def answer(invocation: InvocationOnMock): java.util.Map[String, AnyRef] = {
        delayAndOrThrow()
        new java.util.HashMap[String, AnyRef]()
      }
    }

    def answerFileAttributeView: Answer[FileAttributeView] = new FileSystemAnswer[FileAttributeView](delay, throws) {
      override def answer(invocation: InvocationOnMock): FileAttributeView = {
        delayAndOrThrow()
        mock(classOf[FileAttributeView])
      }
    }

    when(provider.move(any[Path], any[Path])).thenAnswer(answerUnit)
    when(provider.checkAccess(any[Path])).thenAnswer(answerUnit)
    when(provider.createDirectory(any[Path])).thenAnswer(answerUnit)
    when(provider.newByteChannel(any[Path], any[java.util.Set[OpenOption]])).thenAnswer(answerSeekableByteChannel)
    when(provider.isHidden(any[Path])).thenAnswer(answerBoolean)
    when(provider.copy(any[Path], any[Path])).thenAnswer(answerUnit)
    when(provider.delete(any[Path])).thenAnswer(answerUnit)
    when(provider.newDirectoryStream(any[Path], any[Filter[Path]]())).thenAnswer(answerDirectoryStream)
    when(provider.setAttribute(any[Path], any[String], any[Object])).thenAnswer(answerUnit)
    when(provider.readAttributes(any[Path], any[String])).thenAnswer(answerMap)
    when(provider.readAttributes(any[Path], any[Class[BasicFileAttributes]])).thenAnswer(answerBasicFileAttributes)
    when(provider.isSameFile(any[Path], any[Path])).thenAnswer(answerBoolean)
    when(provider.getFileAttributeView(any[Path], any[Class[FileAttributeView]])).thenAnswer(answerFileAttributeView)

    provider
  }

  val testRetryParams = CustomRetryParams.Default.copy(backoff = new Backoff {
    override def next: Backoff = this
    override def backoffMillis: Long = 0
  })

  val pathMock = mock(classOf[Path])

  it should "timeout if the operation takes too long" ignore {
    val retryParams = testRetryParams.copy(timeout = 100 millis)
    val mockFs = mockFileSystem(delay = Option(200 millis))
    val retryableFs = new RetryableFileSystemProviderProxy(mockFs, retryParams)(system)

    a[TimeoutException] shouldBe thrownBy(retryableFs.move(pathMock, pathMock))
    a[TimeoutException] shouldBe thrownBy(retryableFs.checkAccess(pathMock))
    a[TimeoutException] shouldBe thrownBy(retryableFs.createDirectory(pathMock))
    a[TimeoutException] shouldBe thrownBy(retryableFs.newByteChannel(pathMock, mock(classOf[java.util.Set[StandardOpenOption]])))
    a[TimeoutException] shouldBe thrownBy(retryableFs.isHidden(pathMock))
    a[TimeoutException] shouldBe thrownBy(retryableFs.copy(pathMock, pathMock))
    a[TimeoutException] shouldBe thrownBy(retryableFs.delete(pathMock))
    a[TimeoutException] shouldBe thrownBy(retryableFs.newDirectoryStream(pathMock, mock(classOf[Filter[Path]])))
    a[TimeoutException] shouldBe thrownBy(retryableFs.setAttribute(pathMock, "", ""))
    a[TimeoutException] shouldBe thrownBy(retryableFs.readAttributes(pathMock, classOf[BasicFileAttributes]))
    a[TimeoutException] shouldBe thrownBy(retryableFs.readAttributes(pathMock, ""))
    a[TimeoutException] shouldBe thrownBy(retryableFs.isSameFile(pathMock, pathMock))
    a[TimeoutException] shouldBe thrownBy(retryableFs.getFileAttributeView(pathMock, classOf[FileAttributeView]))
  }

  it should "retry on failure and finally succeed if under retry max" in {
    val retryParams = testRetryParams.copy(maxRetries = Option(4))
    val mockFs = mockFileSystem(throws = Option(ThrowParams(new Exception(), nbTimes = 2)))
    val retryableFs = new RetryableFileSystemProviderProxy(mockFs, retryParams)(system)

    retryableFs.move(pathMock, pathMock)
    retryableFs.checkAccess(pathMock)
    retryableFs.createDirectory(pathMock)
    retryableFs.newByteChannel(pathMock, mock(classOf[java.util.Set[StandardOpenOption]]))
    retryableFs.isHidden(pathMock)
    retryableFs.copy(pathMock, pathMock)
    retryableFs.delete(pathMock)
    retryableFs.newDirectoryStream(pathMock, mock(classOf[Filter[Path]]))
    retryableFs.setAttribute(pathMock, "", "")
    retryableFs.readAttributes(pathMock, classOf[BasicFileAttributes])
    retryableFs.readAttributes(pathMock, "")
    retryableFs.isSameFile(pathMock, pathMock)
    retryableFs.getFileAttributeView(pathMock, classOf[FileAttributeView])

    verify(mockFs, times(3)).move(any[Path], any[Path])
    verify(mockFs, times(3)).checkAccess(any[Path])
    verify(mockFs, times(3)).createDirectory(any[Path])
    verify(mockFs, times(3)).newByteChannel(any[Path], any[java.util.Set[OpenOption]])
    verify(mockFs, times(3)).isHidden(any[Path])
    verify(mockFs, times(3)).copy(any[Path], any[Path])
    verify(mockFs, times(3)).delete(any[Path])
    verify(mockFs, times(3)).newDirectoryStream(any[Path], any[Filter[Path]]())
    verify(mockFs, times(3)).setAttribute(any[Path], any[String], any[Object])
    verify(mockFs, times(3)).readAttributes(any[Path], any[String])
    verify(mockFs, times(3)).readAttributes(any[Path], any[Class[BasicFileAttributes]])
    verify(mockFs, times(3)).isSameFile(any[Path], any[Path])
    verify(mockFs, times(3)).getFileAttributeView(any[Path], any[Class[FileAttributeView]])
  }

  it should "retry on failure and fail if over retry max" in {
    val retryParams = testRetryParams.copy(maxRetries = Option(2))
    val mockFs = mockFileSystem(throws = Option(ThrowParams(new IllegalArgumentException(), nbTimes = 3)))
    val retryableFs = new RetryableFileSystemProviderProxy(mockFs, retryParams)(system)

    (the [CromwellFatalException] thrownBy retryableFs.move(pathMock, pathMock)).getCause shouldBe a[IllegalArgumentException]
    (the [CromwellFatalException] thrownBy retryableFs.checkAccess(pathMock)).getCause shouldBe a[IllegalArgumentException]
    (the [CromwellFatalException] thrownBy retryableFs.createDirectory(pathMock)).getCause shouldBe a[IllegalArgumentException]
    (the [CromwellFatalException] thrownBy retryableFs.newByteChannel(pathMock, mock(classOf[java.util.Set[StandardOpenOption]]))).getCause shouldBe a[IllegalArgumentException]
    (the [CromwellFatalException] thrownBy retryableFs.isHidden(pathMock)).getCause shouldBe a[IllegalArgumentException]
    (the [CromwellFatalException] thrownBy retryableFs.copy(pathMock, pathMock)).getCause shouldBe a[IllegalArgumentException]
    (the [CromwellFatalException] thrownBy retryableFs.delete(pathMock)).getCause shouldBe a[IllegalArgumentException]
    (the [CromwellFatalException] thrownBy retryableFs.newDirectoryStream(pathMock, mock(classOf[Filter[Path]]))).getCause shouldBe a[IllegalArgumentException]
    (the [CromwellFatalException] thrownBy retryableFs.setAttribute(pathMock, "", "")).getCause shouldBe a[IllegalArgumentException]
    (the [CromwellFatalException] thrownBy retryableFs.readAttributes(pathMock, classOf[BasicFileAttributes])).getCause shouldBe a[IllegalArgumentException]
    (the [CromwellFatalException] thrownBy retryableFs.readAttributes(pathMock, "")).getCause shouldBe a[IllegalArgumentException]
    (the [CromwellFatalException] thrownBy retryableFs.isSameFile(pathMock, pathMock)).getCause shouldBe a[IllegalArgumentException]
    (the [CromwellFatalException] thrownBy retryableFs.getFileAttributeView(pathMock, classOf[FileAttributeView])).getCause shouldBe a[IllegalArgumentException]

    verify(mockFs, times(3)).move(any[Path], any[Path])
    verify(mockFs, times(3)).checkAccess(any[Path])
    verify(mockFs, times(3)).createDirectory(any[Path])
    verify(mockFs, times(3)).newByteChannel(any[Path], any[java.util.Set[OpenOption]])
    verify(mockFs, times(3)).isHidden(any[Path])
    verify(mockFs, times(3)).copy(any[Path], any[Path])
    verify(mockFs, times(3)).delete(any[Path])
    verify(mockFs, times(3)).newDirectoryStream(any[Path], any[Filter[Path]]())
    verify(mockFs, times(3)).setAttribute(any[Path], any[String], any[Object])
    verify(mockFs, times(3)).readAttributes(any[Path], any[String])
    verify(mockFs, times(3)).readAttributes(any[Path], any[Class[BasicFileAttributes]])
    verify(mockFs, times(3)).isSameFile(any[Path], any[Path])
    verify(mockFs, times(3)).getFileAttributeView(any[Path], any[Class[FileAttributeView]])
  }

  it should "ignore transient exceptions" in {
    def isTransient(t: Throwable) = t.isInstanceOf[FileNotFoundException]
    val retryParams = testRetryParams.copy(maxRetries = Option(1), isTransient = isTransient)
    val mockFs = mockFileSystem(throws = Option(ThrowParams(new FileNotFoundException(), nbTimes = 2)))
    val retryableFs = new RetryableFileSystemProviderProxy(mockFs, retryParams)(system)

    retryableFs.move(pathMock, pathMock)
    retryableFs.checkAccess(pathMock)
    retryableFs.createDirectory(pathMock)
    retryableFs.newByteChannel(pathMock, mock(classOf[java.util.Set[StandardOpenOption]]))
    retryableFs.isHidden(pathMock)
    retryableFs.copy(pathMock, pathMock)
    retryableFs.delete(pathMock)
    retryableFs.newDirectoryStream(pathMock, mock(classOf[Filter[Path]]))
    retryableFs.setAttribute(pathMock, "", "")
    retryableFs.readAttributes(pathMock, classOf[BasicFileAttributes])
    retryableFs.readAttributes(pathMock, "")
    retryableFs.isSameFile(pathMock, pathMock)
    retryableFs.getFileAttributeView(pathMock, classOf[FileAttributeView])

    verify(mockFs, times(3)).move(any[Path], any[Path])
    verify(mockFs, times(3)).checkAccess(any[Path])
    verify(mockFs, times(3)).createDirectory(any[Path])
    verify(mockFs, times(3)).newByteChannel(any[Path], any[java.util.Set[OpenOption]])
    verify(mockFs, times(3)).isHidden(any[Path])
    verify(mockFs, times(3)).copy(any[Path], any[Path])
    verify(mockFs, times(3)).delete(any[Path])
    verify(mockFs, times(3)).newDirectoryStream(any[Path], any[Filter[Path]]())
    verify(mockFs, times(3)).setAttribute(any[Path], any[String], any[Object])
    verify(mockFs, times(3)).readAttributes(any[Path], any[String])
    verify(mockFs, times(3)).readAttributes(any[Path], any[Class[BasicFileAttributes]])
    verify(mockFs, times(3)).isSameFile(any[Path], any[Path])
    verify(mockFs, times(3)).getFileAttributeView(any[Path], any[Class[FileAttributeView]])
  }

  it should "fail immediately on fatal exceptions" in {
    def isFatal(t: Throwable) = t.isInstanceOf[FileNotFoundException]
    val retryParams = testRetryParams.copy(maxRetries = Option(5), isFatal = isFatal)
    val mockFs = mockFileSystem(throws = Option(ThrowParams(new FileNotFoundException(), nbTimes = 3)))
    val retryableFs = new RetryableFileSystemProviderProxy(mockFs, retryParams)(system)

    (the [CromwellFatalException] thrownBy retryableFs.move(pathMock, pathMock)).getCause shouldBe a[FileNotFoundException]
    (the [CromwellFatalException] thrownBy retryableFs.checkAccess(pathMock)).getCause shouldBe a[FileNotFoundException]
    (the [CromwellFatalException] thrownBy retryableFs.createDirectory(pathMock)).getCause shouldBe a[FileNotFoundException]
    (the [CromwellFatalException] thrownBy retryableFs.newByteChannel(pathMock, mock(classOf[java.util.Set[StandardOpenOption]]))).getCause shouldBe a[FileNotFoundException]
    (the [CromwellFatalException] thrownBy retryableFs.isHidden(pathMock)).getCause shouldBe a[FileNotFoundException]
    (the [CromwellFatalException] thrownBy retryableFs.copy(pathMock, pathMock)).getCause shouldBe a[FileNotFoundException]
    (the [CromwellFatalException] thrownBy retryableFs.delete(pathMock)).getCause shouldBe a[FileNotFoundException]
    (the [CromwellFatalException] thrownBy retryableFs.newDirectoryStream(pathMock, mock(classOf[Filter[Path]]))).getCause shouldBe a[FileNotFoundException]
    (the [CromwellFatalException] thrownBy retryableFs.setAttribute(pathMock, "", "")).getCause shouldBe a[FileNotFoundException]
    (the [CromwellFatalException] thrownBy retryableFs.readAttributes(pathMock, classOf[BasicFileAttributes])).getCause shouldBe a[FileNotFoundException]
    (the [CromwellFatalException] thrownBy retryableFs.readAttributes(pathMock, "")).getCause shouldBe a[FileNotFoundException]
    (the [CromwellFatalException] thrownBy retryableFs.isSameFile(pathMock, pathMock)).getCause shouldBe a[FileNotFoundException]
    (the [CromwellFatalException] thrownBy retryableFs.getFileAttributeView(pathMock, classOf[FileAttributeView])).getCause shouldBe a[FileNotFoundException]

    verify(mockFs, times(1)).move(any[Path], any[Path])
    verify(mockFs, times(1)).checkAccess(any[Path])
    verify(mockFs, times(1)).createDirectory(any[Path])
    verify(mockFs, times(1)).newByteChannel(any[Path], any[java.util.Set[OpenOption]])
    verify(mockFs, times(1)).isHidden(any[Path])
    verify(mockFs, times(1)).copy(any[Path], any[Path])
    verify(mockFs, times(1)).delete(any[Path])
    verify(mockFs, times(1)).newDirectoryStream(any[Path], any[Filter[Path]]())
    verify(mockFs, times(1)).setAttribute(any[Path], any[String], any[Object])
    verify(mockFs, times(1)).readAttributes(any[Path], any[String])
    verify(mockFs, times(1)).readAttributes(any[Path], any[Class[BasicFileAttributes]])
    verify(mockFs, times(1)).isSameFile(any[Path], any[Path])
    verify(mockFs, times(1)).getFileAttributeView(any[Path], any[Class[FileAttributeView]])
  }

}
