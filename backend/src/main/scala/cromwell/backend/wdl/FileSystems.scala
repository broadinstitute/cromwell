package cromwell.backend.wdl

import java.nio.file.{FileSystem, Path}

import cromwell.core.PathFactory

trait FileSystems extends PathFactory {

  /**
    * Ordered list of filesystems to be used to execute wdl functions needing IO.
    */
  def fileSystems: List[FileSystem]

  /**
    * Function applied after a string is successfully resolved to a java.nio.Path
    */
  def postMapping(path: Path): Path = path

  /**
    * Use fileSystems in order to try to create a java.nio.Path from path that will be used to perform IO.
    * If no filesystem is able to construct a Path from the String, an exception will be raised.
    */
  protected final def toPath(path: String) = postMapping(buildPath(path, fileSystems))

}
