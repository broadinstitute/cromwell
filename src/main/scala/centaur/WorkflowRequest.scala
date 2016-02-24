package centaur

import java.nio.file.Path

object WorkflowRequest {
  /**
    * Assumes a path like /foo/bar/blah where there'll be blah.wdl, blah.inputs, blah.options. Not
    * at all bulletproof at the moment and will barf all over your face if you do it wrong
    *
    * TODO: It'd be good to make inputs and options optional, I'm pretty sure the submission endpoint views them
    * as optional anyways
    */
  def apply(path: Path): WorkflowRequest = {
    val name = path.getFileName

    val base = path.resolve(name)
    val wdl = base.slurpExtension("wdl")
    val inputs = base.slurpExtension("inputs")
    val options = base.slurpExtension("options")

    WorkflowRequest(name.toString, wdl, inputs, options)
  }

  implicit class EnhancedPath(val path: Path) extends AnyVal {
    def addExtension(extension: String): Path = path.resolveSibling(s"${path.getFileName}.$extension")
    def slurpExtension(extension: String): String = path.addExtension(extension).slurp
    /** Read an entire file into a string, closing the underlying stream. */
    def slurp: String = {
      val source = io.Source.fromFile(path.toFile)
      try source.mkString finally source.close()
    }
  }
}

case class WorkflowRequest(name: String, wdl: String, inputs: String, options: String)
