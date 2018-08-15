package cromwell.core.path

import java.io.InputStream
import java.nio.file.{FileAlreadyExistsException, Files}
import java.nio.file.attribute.{PosixFilePermission, PosixFilePermissions}

import better.files.File.OpenOptions

import scala.collection.JavaConverters._
import scala.concurrent.ExecutionContext
import scala.io.Codec

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

  def createPermissionedDirectories(): this.type = {
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
}
