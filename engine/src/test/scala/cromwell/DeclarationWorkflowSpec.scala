package cromwell

import wdl.{FullyQualifiedName, ImportResolver, WdlNamespaceWithWorkflow}
import cromwell.util.SampleWdl
import org.scalatest.{Matchers, WordSpecLike}
import wom.callable.Callable.{InputDefinitionWithDefault, OptionalInputDefinition, RequiredInputDefinition}
import wom.types.{WomFileType, WomOptionalType, WomStringType}

class DeclarationWorkflowSpec extends Matchers with WordSpecLike {
  "A workflow with declarations in it" should {
    "compute inputs properly" in {

      val expectedRequiredInputs = Map(
        "two_step.cat.file" -> RequiredInputDefinition("two_step.cat.file", WomFileType),
        "two_step.cgrep.str_decl" -> RequiredInputDefinition("two_step.cgrep.str_decl", WomStringType),
        "two_step.cgrep.pattern" -> RequiredInputDefinition("two_step.cgrep.pattern", WomStringType),
        "two_step.flags_suffix" -> RequiredInputDefinition("two_step.flags_suffix", WomStringType)
      )

      /*
       * WARNING: be aware that `workflow.inputs` is used by projects external to Cromwell (eg FC's input enumerator).
       */
      val actualInputs = WdlNamespaceWithWorkflow.load(SampleWdl.DeclarationsWorkflow.workflowSource(), Seq.empty[ImportResolver]).get.workflow.inputs

      actualInputs foreach {
        case (inputName: FullyQualifiedName, actualRequiredInputDefinition: RequiredInputDefinition) =>
          expectedRequiredInputs.get(inputName) match {
            case Some(expectedInputDefinition) => expectedInputDefinition should be(actualRequiredInputDefinition)
            case None => fail(s"Unexpected actual input definition was generated: $inputName")
          }
        case (inputName: FullyQualifiedName, inputDefinition: OptionalInputDefinition) =>
          inputDefinition.localName.value should be("two_step.cat.flags2")
          inputName should be("two_step.cat.flags2")
          inputDefinition.womType should be(WomOptionalType(WomStringType))
        case (inputName: FullyQualifiedName, inputDefinition: InputDefinitionWithDefault) =>
          inputDefinition.localName.value should be ("two_step.static_string")
          inputName should be("two_step.static_string")
          inputDefinition.womType should be(WomStringType)
      }
    }
  }
}
