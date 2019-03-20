package cromwell.services.womtool.models

import io.circe.{Decoder, Encoder, HCursor}
import io.circe.generic.semiauto.deriveEncoder
import wom.callable.Callable.{InputDefinition, InputDefinitionWithDefault, OutputDefinition}
import wom.callable.{CallableTaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle

// Very provisional types for some of these, and perhaps the defaults will go away later in development
case class WorkflowDescription(
                                valid: Boolean,
                                errors: List[String],
                                name: String = "",
                                inputs: List[InputDescription] = List.empty,
                                outputs: List[OutputDescription] = List.empty,
                                images: List[String] = List.empty,
                                submittedDescriptorType: Map[String, String] = Map.empty,
                                importedDescriptorTypes: List[Map[String, String]] = List.empty,
                                meta: Map[String, String] = Map.empty
                              )

case object WorkflowDescription {

  def withErrors(errors: List[String]): WorkflowDescription = {
    WorkflowDescription(
      valid = false,
      errors = errors
    )
  }

  def fromBundle(bundle: WomBundle, languageName: String, languageVersionName: String): WorkflowDescription = {

    val sdt = Map(
      "descriptorType" -> languageName,
      "descriptorTypeVersion" -> languageVersionName
    )

    def createDescription(name: String, inputs: List[InputDefinition], outputs: List[OutputDefinition], meta: Map[String, String]): WorkflowDescription = {
      val inputDescriptions = inputs.sortBy(_.name) map { input: InputDefinition =>
        input match {
          case i: InputDefinitionWithDefault =>
            InputDescription(
              input.name,
              input.womType,
              input.womType.friendlyName,
              input.optional,
              Option(i.default)
            )
          case _ =>
            InputDescription(
              input.name,
              input.womType,
              input.womType.friendlyName,
              input.optional,
              None
            )
        }

      }

      val outputDescriptions = outputs.sortBy(_.name) map { output =>
        OutputDescription(
          output.name,
          output.womType,
          output.womType.friendlyName
        )
      }

      WorkflowDescription(
        valid = true,
        errors = List.empty,
        name = name,
        inputs = inputDescriptions,
        outputs = outputDescriptions,
        images = List.empty,
        submittedDescriptorType = sdt,
        importedDescriptorTypes = List.empty,
        meta = meta
      )
    }

    // Does a source file containing a single task get a primaryCallable? In WDL 1.0 yes, in draft-2 no.
    // https://github.com/broadinstitute/cromwell/pull/3772

    (bundle.allCallables.values.toList, bundle.primaryCallable) match {

      // We have present a workflow or task that this language considers a primary callable
      case (_, Some(primaryCallable: WorkflowDefinition)) =>
        createDescription(primaryCallable.name, primaryCallable.inputs, primaryCallable.outputs, primaryCallable.meta)

      // We have present a workflow or task that this language considers a primary callable
      case (_, Some(primaryCallable: CallableTaskDefinition)) =>
        createDescription(primaryCallable.name, primaryCallable.inputs, primaryCallable.outputs, primaryCallable.meta)

      // WDL draft-2: a solo task is not primary, but we should still use its name and IO
      case ((soloNonPrimaryTask: CallableTaskDefinition) :: Nil, None) =>
        createDescription(soloNonPrimaryTask.name, soloNonPrimaryTask.inputs, soloNonPrimaryTask.outputs, soloNonPrimaryTask.meta)

      // Multiple tasks
      case _ =>
        WorkflowDescription(
          valid = true,
          errors = List.empty,
          name = "", // No name if multiple tasks
          inputs = List.empty,
          outputs = List.empty,
          images = List.empty,
          submittedDescriptorType = sdt,
          importedDescriptorTypes = List.empty,
          meta = Map.empty
        )

    }
  }

  implicit val workflowDescriptionEncoder: Encoder[WorkflowDescription] = deriveEncoder[WorkflowDescription]

  // We need this decoder to exist for `responseAs[WorkflowDescription]` to work in `cromwell.webservice.routes.WomtoolRouteSupportSpec`
  // That test only inspects some fields in the JSON, so this works adequately for now.
  implicit val workflowDescriptionDecoder: Decoder[WorkflowDescription] = (c: HCursor) => {
    for {
      valid <- c.downField("valid").as[Boolean]
      errors <- c.downField("errors").as[List[String]]
    } yield {
      WorkflowDescription(
        valid = valid,
        errors = errors
      )
    }
  }
}
