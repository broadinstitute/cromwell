package wdl4s.cwl

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

case class WorkflowStepId private (fileName: String, stepId: String) extends FullyQualifiedName

object WorkflowStepId {
  def apply(in: String): WorkflowStepId = {
    val Array(fileName, id) = in.split("#")

    WorkflowStepId(fileName, id)
  }
}

case class WorkflowStepOutputId private (fileName: String, stepId: String, outputId: String) extends FullyQualifiedName

object WorkflowStepOutputId {
  def apply(in: String): WorkflowStepOutputId = {
    val Array(fileName, id) = in.split("#")

    val Array(stepId, outputId) = id.split("/")

    WorkflowStepOutputId(fileName, stepId, outputId)
  }
}


case class WorkflowStepOutputIdReference(fileName: String, stepOutputId: String, stepId: String) extends FullyQualifiedName

object WorkflowStepOutputIdReference {
  def apply(in: String): WorkflowStepOutputIdReference = {
    val Array(fileName, stepAndid) = in.split("#")
    val Array(step, id) = stepAndid.split("/")

    WorkflowStepOutputIdReference(fileName, id, step)
  }
}

sealed trait RunOutputId {
  def fileName: String
  def outputId: String
}

case class SameFileRunOutputId(fileName: String, outputId: String, uuid: String, stepId: String) extends RunOutputId

case class DifferentFileRunOutputId(fileName: String, outputId: String) extends RunOutputId

object RunOutputId {
  def apply(in: String): RunOutputId = {
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
     WorkflowStepOutputIdReference(in)
    else
     WorkflowInputId(in)
  }
}
