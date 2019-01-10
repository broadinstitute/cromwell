package cromwell.services.womtool.models

import wom.executable.WomBundle

// Very provisional types for some of these, and perhaps the defaults will go away later in development
case class WorkflowDescription(
                                valid: Boolean,
                                errors: List[String],
                                name: String = "",
                                inputs: List[String] = List.empty,
                                outputs: List[String] = List.empty,
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
        val inputs = callable.inputs map { input =>
          input.name + ", " + input.womType + ", " + (if (input.optional) "optional" else "required")
        }

        val outputs = callable.outputs map { output =>
          output.name + ", " + output.womType
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
          name = "",
          inputs = List.empty,
          outputs = List.empty,
          images = List.empty,
          submittedDescriptorType = sdt,
          importedDescriptorTypes = List.empty,
          meta = Map.empty
        )

    }
  }

}
