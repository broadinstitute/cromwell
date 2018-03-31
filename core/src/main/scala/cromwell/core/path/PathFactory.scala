package cromwell.core.path

import cromwell.core.path.PathFactory.PathBuilders

import scala.util.Success

/**
  * Convenience trait delegating to the PathFactory singleton
  */
trait PathFactory {
  /**
    * Path builders to be applied (in order) to attempt to build a Path from a string.
    */
  def pathBuilders: PathBuilders

  /**
    * Function applied after a string is successfully resolved to a Path
    */
  def postMapping(path: Path): Path = path

  /**
    * Function applied before a string is attempted to be resolved to a Path
    */
  def preMapping(string: String): String = string

  /**
    * Attempts to build a Path from a String
    */
  def buildPath(string: String): Path = PathFactory.buildPath(string, pathBuilders, preMapping, postMapping)
}

object PathFactory {
  type PathBuilders = List[PathBuilder]

  /**
    * Attempts to build a Path from a String
    */
  def buildPath(string: String,
                pathBuilders: List[PathBuilder],
                preMapping: String => String = identity[String],
                postMapping: Path => Path = identity[Path]): Path = {
    pathBuilders.toStream map { _.build(preMapping(string)) } collectFirst { case Success(p) => postMapping(p) } getOrElse {
      val pathBuilderNames: String = pathBuilders map { _.name } mkString ", "
      throw PathParsingException(
        s"$string exists on a filesystem not supported by this instance of Cromwell." +
        s" Supported filesystems are: $pathBuilderNames." +
        s" Please refer to the documentation for more information on how to configure filesystems: http://cromwell.readthedocs.io/en/develop/backends/HPC/#filesystems"
      )
    }
  }
}
