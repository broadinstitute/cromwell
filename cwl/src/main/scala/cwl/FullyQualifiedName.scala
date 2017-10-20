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
}

case class WorkflowInputId private (fileName: String, inputId: String) extends FullyQualifiedName

object WorkflowInputId {
  def apply(in: String): WorkflowInputId = {
    val Array(fileName, id) = in.split("#")

    WorkflowInputId(fileName, id)
  }
}

case class WorkflowOutputId private (fileName: String, stepId: String, outputId: String) extends FullyQualifiedName

object WorkflowOutputId {
  def apply(in: String): WorkflowOutputId = {
    val Array(fileName, id) = in.split("#")
    val Array(stepId, outputId) = id.split("/")

    WorkflowOutputId(fileName, stepId, outputId)
  }
}

case class WorkflowStepId private (fileName: String, stepId: String) extends FullyQualifiedName

object WorkflowStepId {
  def apply(in: String): WorkflowStepId = {
    val Array(fileName, id) = in.split("#")

    WorkflowStepId(fileName, id)
  }
}

case class WorkflowStepInputOrOutputId private(fileName: String, stepId: String, ioId: String) extends FullyQualifiedName

object WorkflowStepInputOrOutputId {
  def apply(in: String): WorkflowStepInputOrOutputId = {
    val Array(fileName, id) = in.split("#")
    val Array(stepId, fieldId) = id.split("/")

    WorkflowStepInputOrOutputId(fileName, stepId, fieldId)
  }
}

sealed trait RunId {
  def fileName: String
  def variableId: String
}

case class SameFileRunOutputId(fileName: String, variableId: String, uuid: String, stepId: String) extends RunId

case class DifferentFileRunOutputId(fileName: String, variableId: String) extends RunId

object RunId {
  def apply(in: String): RunId = {
    val Array(fileName, stepAndid) = in.split("#")

    if (stepAndid.contains("/")) {
      val Array(step, uuid, id) = stepAndid.split("/")
      SameFileRunOutputId(fileName, id, uuid, step)
    } else
      DifferentFileRunOutputId(fileName, stepAndid)
  }
}
object FullyQualifiedName {
  def apply(in: String): FullyQualifiedName = {

   val Array(_, after) = in.split("#")

   if (after.contains("/"))
     WorkflowStepInputOrOutputId(in)
    else
     WorkflowInputId(in)
  }
}
