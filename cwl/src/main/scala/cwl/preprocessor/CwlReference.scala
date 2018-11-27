package cwl.preprocessor

import better.files.{File => BFile}
import cwl.preprocessor.CwlReference._
import cwl.{FileAndId, FullyQualifiedName}
import cwl.command.ParentName

sealed trait CwlReference {
  def pathAsString: String
  def pointerWithinFile: Option[String]
  def changePointer(to: Option[String]): CwlReference

  private def pointerWithHash: String = pointerWithinFile.map(p => s"#$p").getOrElse("")
  lazy val fullReference: String = s"$pathAsString$pointerWithHash"

  override def toString: String = fullReference
}

/**
  * Saladed CWLs reference other local CWL "node" (workflow or tool) using a URI as follow:
  * file:///path/to/file/containing/node.cwl[#pointer_to_node]
  * #pointer_to_node to node is optional, and will specify which workflow or tool is being targeted in the file.
  *
  * e.g:
  *   {
  *     "class": "Workflow",
  *     "id": "file:///path/to/workflow/workflow.cwl",
  *     ...
  *     "steps": [
  *       {
  *         "run": "file:///path/to/workflow/multi_tools.cwl#my_tool",
  *         ...
  *       }
  *     ]
  *   }
  *
  * This snippet contains 2 references, one that is the ID of this workflow, the other one is the run step pointing to "my_tool" in "/path/to/workflow/multi_tools.cwl"
  *
  */
final case class CwlFileReference(file: BFile, pointerWithinFile: Option[String]) extends CwlReference {
  override val pathAsString: String = s"$LocalScheme${file.toString}"
  override def changePointer(to: Option[String]): CwlReference = this.copy(pointerWithinFile = to)
}

final case class CwlHttpReference(pathAsString: String, pointerWithinFile: Option[String]) extends CwlReference {
  override def changePointer(to: Option[String]): CwlReference = this.copy(pointerWithinFile = to)
}

object CwlReference {
  val LocalScheme = "file://"
  val HttpScheme = "http://"
  val HttpsScheme = "https://"

  implicit class EnhancedCwlId(val id: String) extends AnyVal {
    def asReference: Option[CwlReference] = CwlReference.fromString(id)
    def stripFilePrefix = id.stripPrefix(LocalScheme)
  }

  val ReferenceRegex = "(.*://)?([^#]*)(#(.*))?".r

  def fromString(in: String): Option[CwlReference] = {
    in match {
      case ReferenceRegex(scheme, path, _, pointerWithinFile) =>
        if (scheme == LocalScheme) {
          FullyQualifiedName.maybeApply(in)(ParentName.empty) match {
            case Some(FileAndId(file, _, _)) => Option(CwlFileReference(BFile(file.stripFilePrefix), Option(pointerWithinFile)))
            case _ => Option(CwlFileReference(BFile(in.stripFilePrefix), Option(pointerWithinFile)))
          }
        } else if (scheme == HttpScheme || scheme == HttpsScheme) {
          Option(CwlHttpReference(s"$scheme$path", Option(pointerWithinFile)))
        } else {
          None
        }
    }
  }
}

object CwlFileReference {
  def apply(file: BFile, pointer: Option[String]) = {
    new CwlFileReference(file, pointer)
  }
}
