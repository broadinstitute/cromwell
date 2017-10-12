package cwl

/**
  * All of these classes decompose a "fully qualified" id into its constitutent parts.  They have unique types as
  * they are seen in different parts of a CWL document and they differ in content.
  *
  * The fully qualified names are created by the Schema salad preprocessing step.
  *
  * @see <a href="http://www.commonwl.org/v1.0/SchemaSalad.html#Identifier_resolution">Schema salad Identifier Resolution</a>
  */
trait FullyQualifiedName {
  val fileName: String
  val id: String
}

case class FileAndId private(fileName: String, id: String) extends FullyQualifiedName

object FileAndId {
  def apply(in: String): FileAndId = {
    val Array(fileName, id) = in.split("#")

    FileAndId(fileName, id)
  }
}

case class FileStepAndId private(fileName: String, stepId: String, id: String) extends FullyQualifiedName

object FileStepAndId {
  def apply(in: String): FileStepAndId = {
    val Array(fileName, id) = in.split("#")
    val Array(stepId, outputId) = id.split("/")

    FileStepAndId(fileName, stepId, outputId)
  }
}

case class FileStepUUID(fileName: String, id: String, uuid: String, stepId: String) extends FullyQualifiedName

object FullyQualifiedName {
  def apply(in: String): FullyQualifiedName = {

    val Array(file, after) = in.split("#")

    (after.split("/").toList) match {
      case step :: uuid :: id :: Nil => FileStepUUID(file, id, uuid, step)
      case step :: id :: Nil => FileStepAndId(file, step, id)
      case id :: Nil => FileAndId(file, id)
      case _ => throw new RuntimeException(s"malformed FQN: $in")
    }
  }
}
