package cromwell.services.womtool.models

import cromwell.services.womtool.models.InputDescription.inputDescriptionWriter
import cromwell.services.womtool.models.MetaValueElementJsonSupport.metaValueElementEncoder
import cromwell.services.womtool.models.OutputDescription.outputDescriptionWriter
import spray.json.DefaultJsonProtocol._
import spray.json._
import wom.RuntimeAttributesKeys
import wom.callable.Callable.{
  InputDefinition,
  InputDefinitionWithDefault,
  OutputDefinition,
  RuntimeOverrideInputDefinition
}
import wom.callable.{Callable, CallableTaskDefinition, MetaValueElement, WorkflowDefinition}
import wom.executable.WomBundle

case class WorkflowDescription(valid: Boolean,
                               errors: List[String],
                               validWorkflow: Boolean,
                               name: String,
                               inputs: List[InputDescription],
                               outputs: List[OutputDescription],
                               images: List[String],
                               submittedDescriptorType: Map[String, String],
                               importedDescriptorTypes: List[Map[String, String]],
                               meta: Map[String, MetaValueElement],
                               parameterMeta: Map[String, MetaValueElement],
                               isRunnableWorkflow: Boolean
)

case object WorkflowDescription {

  def withErrors(errors: List[String], submittedDescriptorType: Map[String, String]): WorkflowDescription =
    WorkflowDescription(
      valid = false,
      errors = errors,
      validWorkflow = false,
      submittedDescriptorType = submittedDescriptorType
    )

  def fromBundle(bundle: WomBundle,
                 submittedDescriptorType: Map[String, String],
                 inputErrors: List[String] = List.empty
  ): WorkflowDescription = {

    val images: List[String] = bundle.allCallables.values.toList flatMap { callable: Callable =>
      callable match {
        case task: CallableTaskDefinition =>
          task.runtimeAttributes.attributes.get(RuntimeAttributesKeys.DockerKey).map(_.sourceString)
        case _ =>
          None
      }
    }

    // Does a source file containing a single task get a primaryCallable? In WDL 1.0 yes, in draft-2 no.
    // https://github.com/broadinstitute/cromwell/pull/3772

    (bundle.allCallables.values.toList, bundle.primaryCallable) match {

      // There is a primary callable in the form of a workflow
      case (_, Some(primaryCallable: WorkflowDefinition)) =>
        fromBundleInner(
          inputErrors,
          primaryCallable.name,
          submittedDescriptorType,
          primaryCallable.inputs,
          primaryCallable.outputs,
          primaryCallable.meta,
          primaryCallable.parameterMeta,
          images,
          isRunnableWorkflow = true
        )

      // There is a primary callable in the form of a task
      case (_, Some(primaryCallable: CallableTaskDefinition)) =>
        fromBundleInner(
          inputErrors,
          primaryCallable.name,
          submittedDescriptorType,
          primaryCallable.inputs,
          primaryCallable.outputs,
          primaryCallable.meta,
          primaryCallable.parameterMeta,
          images,
          isRunnableWorkflow = false
        )

      // WDL draft-2: a solo task is not primary, but we should still use its name and IO
      case ((soloNonPrimaryTask: CallableTaskDefinition) :: Nil, None) =>
        fromBundleInner(
          inputErrors,
          soloNonPrimaryTask.name,
          submittedDescriptorType,
          soloNonPrimaryTask.inputs,
          soloNonPrimaryTask.outputs,
          soloNonPrimaryTask.meta,
          soloNonPrimaryTask.parameterMeta,
          images,
          isRunnableWorkflow = false
        )

      // Multiple tasks
      case _ =>
        fromBundleInner(inputErrors,
                        "",
                        submittedDescriptorType,
                        List.empty,
                        List.empty,
                        Map.empty,
                        Map.empty,
                        images,
                        isRunnableWorkflow = false
        )
    }
  }

  private def fromBundleInner(
    inputErrors: List[String],
    name: String,
    submittedDescriptorType: Map[String, String],
    inputs: List[InputDefinition],
    outputs: List[OutputDefinition],
    meta: Map[String, MetaValueElement],
    parameterMeta: Map[String, MetaValueElement],
    images: List[String],
    isRunnableWorkflow: Boolean
  ): WorkflowDescription = {
    val inputDescriptions = inputs.sortBy(_.name) flatMap { input: InputDefinition =>
      input match {
        case _: RuntimeOverrideInputDefinition =>
          None // Exclude runtime override inputs from the description
        case i: InputDefinitionWithDefault =>
          Option(
            InputDescription(
              input.name,
              input.womType,
              input.womType.friendlyName,
              input.optional,
              Option(i.default)
            )
          )
        case _ =>
          Option(
            InputDescription(
              input.name,
              input.womType,
              input.womType.friendlyName,
              input.optional,
              None
            )
          )
      }

    }

    val outputDescriptions = outputs.sortBy(_.name) map { output: OutputDefinition =>
      OutputDescription(
        output.name,
        output.womType,
        output.womType.friendlyName
      )
    }

    new WorkflowDescription(
      valid = inputErrors.isEmpty,
      errors = inputErrors,
      validWorkflow = true,
      name = name,
      inputs = inputDescriptions,
      outputs = outputDescriptions,
      images = images,
      submittedDescriptorType = submittedDescriptorType,
      importedDescriptorTypes = List.empty,
      meta = meta,
      parameterMeta = parameterMeta,
      isRunnableWorkflow = isRunnableWorkflow
    )
  }

  def apply(valid: Boolean,
            errors: List[String],
            validWorkflow: Boolean,
            name: String = "",
            inputs: List[InputDescription] = List.empty,
            outputs: List[OutputDescription] = List.empty,
            images: List[String] = List.empty,
            submittedDescriptorType: Map[String, String] = Map.empty,
            importedDescriptorTypes: List[Map[String, String]] = List.empty,
            meta: Map[String, MetaValueElement] = Map.empty,
            parameterMeta: Map[String, MetaValueElement] = Map.empty,
            isRunnableWorkflow: Boolean = false
  ): WorkflowDescription =
    new WorkflowDescription(valid,
                            errors,
                            validWorkflow,
                            name,
                            inputs,
                            outputs,
                            images,
                            submittedDescriptorType,
                            importedDescriptorTypes,
                            meta,
                            parameterMeta,
                            isRunnableWorkflow
    )

  implicit val workflowDescriptionFormat: RootJsonFormat[WorkflowDescription] =
    new RootJsonFormat[WorkflowDescription] {
      override def write(a: WorkflowDescription): JsValue = JsObject(
        "valid" -> JsBoolean(a.valid),
        "errors" -> JsArray(a.errors.map(JsString(_)).toVector),
        "validWorkflow" -> JsBoolean(a.validWorkflow),
        "name" -> JsString(a.name),
        "inputs" -> JsArray(a.inputs.map(inputDescriptionWriter.write).toVector),
        "outputs" -> JsArray(a.outputs.map(outputDescriptionWriter.write).toVector),
        "images" -> JsArray(a.images.map(JsString(_)).toVector),
        "submittedDescriptorType" -> stringMap(a.submittedDescriptorType),
        "importedDescriptorTypes" -> JsArray(a.importedDescriptorTypes.map(stringMap).toVector),
        "meta" -> metaMap(a.meta),
        "parameterMeta" -> metaMap(a.parameterMeta),
        "isRunnableWorkflow" -> JsBoolean(a.isRunnableWorkflow)
      )

      // We need this reader to exist for `responseAs[WorkflowDescription]` to work in `cromwell.webservice.routes.WomtoolRouteSupportSpec`
      // That test only inspects some fields in the JSON, so this works adequately for now.
      override def read(json: JsValue): WorkflowDescription = {
        val fields = json.asJsObject.fields
        WorkflowDescription(
          valid = fields("valid").convertTo[Boolean],
          errors = fields("errors").convertTo[List[String]],
          validWorkflow = fields("validWorkflow").convertTo[Boolean]
        )
      }

      private def stringMap(m: Map[String, String]): JsObject = JsObject(m.map { case (k, v) => k -> JsString(v) })
      private def metaMap(m: Map[String, MetaValueElement]): JsObject =
        JsObject(m.map { case (k, v) => k -> metaValueElementEncoder.write(v) })
    }
}
