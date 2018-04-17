package cwl

import cwl.command.ParentName

/**
  * All of these classes decompose a "fully qualified" id into its constituent parts.  They have unique types as
  * they are seen in different parts of a CWL document and they differ in content.
  *
  * The fully qualified names are created by the Schema salad preprocessing step.
  *
  * @see <a href="http://www.commonwl.org/v1.0/SchemaSalad.html#Identifier_resolution">Schema salad Identifier Resolution</a>
  */
trait FullyQualifiedName {
  def fileName: String
  def id: String
  def parent: Option[String]
}

case class FileAndId private(fileName: String, parent: Option[String], id: String) extends FullyQualifiedName

object FileAndId {
  def apply(in: String)(implicit parent: ParentName): FileAndId = {
    val Array(fileName, id) = in.split("#")
    val cleanID = parent.stripParent(id)
    FileAndId(fileName, parent.value, cleanID)
  }
}

case class FileStepAndId private(fileName: String, parent: Option[String], stepId: String, id: String) extends FullyQualifiedName

object FileStepAndId {
  def apply(in: String)(implicit parent: ParentName): FileStepAndId = {
    val Array(fileName, id) = in.split("#")
    val cleanID = parent.stripParent(id)
    val Array(stepId, outputId) = cleanID.split("/")
    FileStepAndId(fileName, parent.value, stepId, outputId)
  }
}

case class ArbitrarilyNested(fileName: String, parent: Option[String], id: String) extends FullyQualifiedName

case class FileStepUUID(fileName: String, parent: Option[String], id: String, uuid: String, stepId: String) extends FullyQualifiedName

object FullyQualifiedName {
  def apply(in: String)(implicit parent: ParentName): FullyQualifiedName = maybeApply(in)(parent).getOrElse(throw new Exception(s"malformed FQN: $in"))

  def maybeApply(in: String)(implicit parent: ParentName): Option[FullyQualifiedName] = {

    in.split("#") match {
      case Array(file, after) =>
        val cleanAfter = parent.stripParent(after)
        cleanAfter.split("/").toList match {
          case step :: uuid :: id :: Nil => Option(FileStepUUID(file, parent.value, id, uuid, step))
          case step :: id :: Nil => Option(FileStepAndId(file, parent.value, step, id))
          case id :: Nil => Option(FileAndId(file, parent.value, id))
          case many => Option(ArbitrarilyNested(file, parent.value, many.last))
        }
      case _ => None
    }
  }
}
