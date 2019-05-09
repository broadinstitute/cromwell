package cromwell.core.path

import java.io.{BufferedReader, IOException, InputStream, InputStreamReader}
import java.nio.file.{FileAlreadyExistsException, Files}
import java.nio.file.attribute.{PosixFilePermission, PosixFilePermissions}

import better.files.File.OpenOptions
import cromwell.util.TryWithResource.tryWithResource

import scala.collection.JavaConverters._
import scala.concurrent.ExecutionContext
import scala.io.Codec
import scala.util.Failure

/**
  * Implements methods beyond those implemented in NioPathMethods and BetterFileMethods
  *
  * @see [[cromwell.core.path.Path]]
  * @see [[cromwell.core.path.NioPathMethods]]
  * @see [[cromwell.core.path.BetterFileMethods]]
  */
trait EvenBetterPathMethods {
  self: Path =>

  final def getAttribute(attribute: String): Any = java.nio.file.Files.getAttribute(nioPathPrivate, attribute)

  final def plusExt(ext: String): Path = plusSuffix(s".$ext")

  final def swapExt(oldExt: String, newExt: String): Path = swapSuffix(s".$oldExt", s".$newExt")

  final def nameWithoutExtensionNoIo: String = if (name contains ".") name.substring(0, name lastIndexOf ".") else name

  final def plusSuffix(suffix: String): Path = swapSuffix("", suffix)

  final def swapSuffix(oldSuffix: String, newSuffix: String): Path = {
    sibling(s"${name.stripSuffix(oldSuffix)}$newSuffix")
  }

  final def createTempFile(prefix: String = "", suffix: String = ""): Path = {
    newPath(java.nio.file.Files.createTempFile(nioPathPrivate, prefix, suffix))
  }

  def chmod(permissions: String): this.type = {
    setPermissions(PosixFilePermissions.fromString(permissions).asScala.toSet)
    this
  }

  // betterFile.symbolicLink calls Files.readSymbolicLink, but then implicitly converts the java.nio.Path returned to a better.File
  // which calls toAbsolutePath. Consequently, if the path was relative, the current directory is used to make it absolute.
  // This is not the desired behavior to be able to follow relative symbolic links, so bypass better files method and directly use the java one.
  final def symbolicLinkRelative: Option[Path] = {
    if (betterFile.isSymbolicLink) {
      Option(newPath(Files.readSymbolicLink(betterFile.path)))
    } else None
  }

  final def followSymbolicLinks: Path = {
    symbolicLinkRelative match {
      case Some(target) => parent.resolve(target.followSymbolicLinks)
      case None => this
    }
  }

  final def createPermissionedDirectories(): this.type = {
    if (!exists) {
      if (parent != null) parent.createPermissionedDirectories()
      try {
        createDirectories()
        // When using PosixFilePermissions/FileAttributes with createDirectories, the umask Cromwell happens to be using
        // affects the resulting directory permissions.  This is not the desired behavior, these directories should be
        // world readable/writable/executable irrespective of the umask.
        addPermission(PosixFilePermission.OTHERS_READ)
        addPermission(PosixFilePermission.OTHERS_WRITE)
        addPermission(PosixFilePermission.OTHERS_EXECUTE)
      }
      catch {
        // Race condition that's particularly likely with scatters.  Ignore.
        case _: FileAlreadyExistsException =>
        // The GCS filesystem does not support setting permissions and will throw an `UnsupportedOperationException`.
        // Evaluating expressions like `write_lines` in a command block will cause the above permission-manipulating
        // code to run against a GCS Path. Fortunately creating directories in GCS is also unnecessary, so this
        // exception type is just ignored.
        case _: UnsupportedOperationException =>
      }
    }
    this
  }

  final def untailed = UntailedWriter(this)

  final def tailed(tailedSize: Int) = TailedWriter(this, tailedSize)

  def mediaInputStream(implicit ec: ExecutionContext): InputStream = {
    // See https://github.com/scala/bug/issues/10347 and https://github.com/scala/bug/issues/10790
    locally(ec)
    newInputStream
  }

  def writeContent(content: String)(openOptions: OpenOptions, codec: Codec)(implicit ec: ExecutionContext): this.type = {
    locally(ec)
    write(content)(openOptions, Codec.UTF8)
  }

  /*
   * The input stream will be closed when this method returns, which means the f function
   * cannot leak an open stream.
   */
  def withReader[A](f: BufferedReader => A)(implicit ec: ExecutionContext): A = {

    // Use an input reader to convert the byte stream to character stream. Buffered reader for efficiency.
    tryWithResource(() => new BufferedReader(new InputStreamReader(this.mediaInputStream, Codec.UTF8.name)))(f).recoverWith({
      case failure => Failure(new IOException(s"Could not read from ${this.pathAsString}: ${failure.getMessage}", failure))
    }).get
  }

  /**
    * Returns an Array[Byte] from a Path. Limit the array size to "limit" byte if defined.
    * @throws IOException if failOnOverflow is true and the file is larger than limit
    */
  def limitFileContent(limit: Option[Int], failOnOverflow: Boolean)(implicit ec: ExecutionContext) = withReader { reader =>
    val bytesIterator = Iterator.continually(reader.read).takeWhile(_ != -1).map(_.toByte)
    // Take 1 more than the limit so that we can look at the size and know if it's overflowing
    val bytesArray = limit.map(l => bytesIterator.take(l + 1)).getOrElse(bytesIterator).toArray

    limit match {
      case Some(l) if failOnOverflow && bytesArray.length > l =>
        throw new IOException(s"File $this is larger than $l Bytes. Maximum read limits can be adjusted in the configuration under system.input-read-limits.")
      case Some(l) => bytesArray.take(l)
      case _ => bytesArray
    }
  }

  /**
    * Reads the first limitBytes of a file and makes a String. Prepend with an annotation at the start (to say that this is the
    * first n bytes).
    */
  def annotatedContentAsStringWithLimit(limitBytes: Int)(implicit ec: ExecutionContext): String =
    s"[First $limitBytes bytes]:" + new String(limitFileContent(Option(limitBytes), failOnOverflow = false))
}
