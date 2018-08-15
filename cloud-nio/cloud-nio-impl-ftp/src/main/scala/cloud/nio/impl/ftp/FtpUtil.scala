package cloud.nio.impl.ftp

import java.io.IOException
import java.nio.file.{FileAlreadyExistsException, NoSuchFileException}

import io.github.andrebeat.pool.Lease
import org.apache.commons.net.ftp.{FTPClient, FTPReply}

import scala.util.{Failure, Success, Try}

object FtpUtil {
  case class FtpIoException(message: String, code: Int, replyString: String, cause: Throwable = null) extends IOException(s"$message: $replyString", cause) {
    def isTransient = FTPReply.isNegativeTransient(code)
    def isFatal = FTPReply.isNegativePermanent(code)
  }
  
  case class FtpOperation(cloudHost: String, cloudPath: String, description: String) {
    def fullPath = s"ftp://$cloudHost/$cloudPath"

    def errorMessage = s"Failed to $description at $fullPath"

    def fail(lease: Lease[FTPClient], cause: Option[Throwable] = None) = {
      val exception = generateException(lease.get(), cause)

      /*
       * The idea is that the ftp client threw an exception, something might be wrong with the connection so invalidate (destroy) it.
       * If the operation failed without throwing, release the client to the pool without trashing it for good.
       */
      if (cause.isDefined) lease.invalidate() else lease.release()
      
      throw exception
    }
    
    def generateException(client: FTPClient, cause: Option[Throwable]) = cause match {
      case None if client.getReplyCode == FTPReply.FILE_UNAVAILABLE && client.getReplyString.toLowerCase.contains("exists") => 
        new FileAlreadyExistsException(fullPath)
      case None if client.getReplyCode == FTPReply.FILE_UNAVAILABLE && client.getReplyString.toLowerCase.contains("no such file") => 
        new NoSuchFileException(fullPath)
      case None => FtpIoException(errorMessage, client.getReplyCode, Option(client.getReplyString).getOrElse("N/A"))
      case Some(c) => FtpIoException(errorMessage, client.getReplyCode, Option(client.getReplyString).getOrElse("N/A"), c)
    }

    override def toString = s"$description at $fullPath"
  }

  def runBoolean(operation: FtpOperation, lease: Lease[FTPClient], autoRelease: Boolean = true)(f: FTPClient => Boolean): Unit = {
    Try(f(lease.get())) match {
      // Operation didn't throw but the result is false which means it failed
      case Success(false) => operation.fail(lease)
      case Success(true) => if (autoRelease) lease.release()
      case Failure(failure) => operation.fail(lease, Option(failure))
    }
  }

  def run[A <: AnyRef](operation: FtpOperation, lease: Lease[FTPClient], autoRelease: Boolean = true)(f: FTPClient => A): A = {
    Try(f(lease.get())) match {
      // Operation didn't throw but the result is null which means it failed
      case Success(null) => operation.fail(lease)
      case Success(result) =>
        if (autoRelease) lease.release()
        result
      case Failure(failure) => operation.fail(lease, Option(failure))
    }
  }
}
