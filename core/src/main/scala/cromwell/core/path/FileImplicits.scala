package cromwell.core.path

import java.nio.file.FileAlreadyExistsException
import java.nio.file.attribute.PosixFilePermission

import better.files.File

object FileImplicits {

  implicit class EnhancedFile(val directory: File) extends AnyVal {
    def createPermissionedDirectories(): directory.type = {
      if (!directory.exists) {
        directory.parent.createPermissionedDirectories()
        try {
          directory.createDirectory()
          // When using PosixFilePermissions/FileAttributes with createDirectories, the umask Cromwell happens to be using
          // affects the resulting directory permissions.  This is not the desired behavior, these directories should be
          // world readable/writable/executable irrespective of the umask.
          directory.addPermission(PosixFilePermission.OTHERS_READ)
          directory.addPermission(PosixFilePermission.OTHERS_WRITE)
          directory.addPermission(PosixFilePermission.OTHERS_EXECUTE)
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
      directory
    }
  }
}
