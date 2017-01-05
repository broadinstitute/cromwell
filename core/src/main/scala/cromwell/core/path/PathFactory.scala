package cromwell.core.path

import java.nio.file.Path

import better.files.File

import scala.util.Success

/**
  * Convenience trait delegating to the PathFactory singleton
  */
trait PathFactory {
  /**
    * Path builders to be applied (in order) to attempt to build a java.nio.Path from a string.
    */
  def pathBuilders: List[PathBuilder]

  /**
    * Function applied after a string is successfully resolved to a java.nio.Path
    */
  def postMapping(path: Path): Path = path

  /**
    * Function applied before a string is attempted to be resolved to a java.nio.Path
    */
  def preMapping(string: String): String = string

  /**
    * Attempts to build a java.nio.Path from a String
    */
  def buildPath(string: String): Path = PathFactory.buildPath(string, pathBuilders, preMapping, postMapping)

  /**
    * Attempts to build a better.files.File from a String
    */
  def buildFile(string: String): File = PathFactory.buildFile(string, pathBuilders, preMapping, postMapping)
}

object PathFactory {
  /**
    * Attempts to build a java.nio.Path from a String
    */
  def buildPath(string: String,
                pathBuilders: List[PathBuilder],
                preMapping: String => String = identity[String],
                postMapping: Path => Path = identity[Path]): Path = {
    pathBuilders.toStream map { _.build(preMapping(string)) } collectFirst { case Success(p) => postMapping(p) } getOrElse {
      val pathBuilderNames: String = pathBuilders map { _.name } mkString ", "
      throw PathParsingException(s"Could not find suitable filesystem among $pathBuilderNames to parse $string.")
    }
  }

  def buildFile(string: String,
                pathBuilders: List[PathBuilder],
                preMapping: String => String = identity[String],
                postMapping: Path => Path = identity[Path]): File = File(buildPath(string, pathBuilders, preMapping, postMapping))

  def pathPlusSuffix(path: File, suffix: String): File = path.sibling(s"${path.name}.$suffix")
}
