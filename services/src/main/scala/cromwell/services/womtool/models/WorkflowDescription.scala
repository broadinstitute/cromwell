package cromwell.services.womtool.models

import io.circe.{Decoder, Encoder, HCursor}
import io.circe.generic.semiauto.deriveEncoder
import wom.callable.Callable.{InputDefinitionWithDefault, InputDefinition}
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

    bundle.primaryCallable match {
      case Some(callable) =>
        val inputs = callable.inputs.sortBy(_.name) map { input: InputDefinition =>
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

        val outputs = callable.outputs.sortBy(_.name) map { output =>
          OutputDescription(
            output.name,
            output.womType,
            output.womType.friendlyName
          )
        }

        WorkflowDescription(
          valid = true,
          errors = List.empty,
          name = callable.name,
          inputs = inputs,
          outputs = outputs,
          images = List.empty,
          submittedDescriptorType = sdt,
          importedDescriptorTypes = List.empty,
          meta = Map.empty
        )

      case None =>
        WorkflowDescription(
          valid = true,
          errors = List.empty,
          name = if (bundle.allCallables.size == 1) bundle.allCallables.head._1 else "", // No name if multiple tasks
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
