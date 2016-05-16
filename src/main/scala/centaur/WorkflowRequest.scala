package centaur

import java.io.FileNotFoundException
import java.nio.file.Path

import centaur.test.standard.StandardTestCase

import scala.util.{Failure, Success, Try}

/**
  * Wraps the necessary information to submit a single workflow request to Cromwell.
  *
  * TODO: base & metadata can be removed when we modify metadata handling to use this new structure
  * TODO: can we factor out name? it's not *really* part of the request, only used for output/exceptions
  */
case class WorkflowRequest(name: String, base: Path, wdl: String, inputs: Option[String], options: Option[String], metadata: Option[String])

object WorkflowRequest {
  // TODO: The slurps will throw - not a high priority but see #36
  def apply(testCase: StandardTestCase): WorkflowRequest = {
    WorkflowRequest(testCase.name,
      testCase.basePath,
      testCase.files.wdl.slurp,
      testCase.files.inputs map { _.slurp },
      testCase.files.options map { _.slurp },
      testCase.files.metadata map { _.slurp }
    )
  }

  /**
    * Assumes a path like /foo/bar/blah where there'll be blah.wdl, blah.inputs, blah.options. Not
    * at all bulletproof at the moment and will barf all over your face if you do it wrong
    */
  def apply(path: Path): WorkflowRequest = {
    val name = path.getFileName
    val base = path.resolve(name)
    val wdl = base.slurpExtension("wdl")
    val inputs = base.slurpExtensionIfExists("inputs")
    val options = base.slurpExtensionIfExists("options")
    val metadata = base.slurpExtensionIfExists("metadata")

    WorkflowRequest(name.toString, base, wdl, inputs, options, metadata)
  }

  implicit class EnhancedPath(val path: Path) extends AnyVal {
    def addExtension(extension: String): Path = path.resolveSibling(s"${path.getFileName}.$extension")
    def slurpExtension(extension: String): String = path.addExtension(extension).slurp

    def slurpExtensionIfExists(extension: String): Option[String] = {
      val attempt = Try(slurpExtension(extension))
      attempt match {
        case Success(x) => Option(x)
        case Failure(t: FileNotFoundException) => None
        case Failure(t) => throw t
      }
    }

    /** Read an entire file into a string, closing the underlying stream. */
    def slurp: String = {
      val source = io.Source.fromFile(path.toFile)
      try source.mkString finally source.close()
    }
  }
}
