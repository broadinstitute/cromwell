package cromwell.api.model

import spray.json.{DefaultJsonProtocol, JsObject, JsValue}

object WorkflowDescriptionJsonSupport extends DefaultJsonProtocol {
  implicit val WaasWorkflowDescriptorTypeFormat = jsonFormat2(WaasWorkflowDescriptorType)
  implicit val WaasDescriptionWomTypeFormat = jsonFormat1(WaasDescriptionWomType)
  implicit val WaasDescriptionInputDefinitionFormat = jsonFormat5(WaasDescriptionInputDefinition)
  implicit val WaasDescriptionOutputDefinitionFormat = jsonFormat3(WaasDescriptionOutputDefinition)

  implicit val WorkflowDescriptionFormat = jsonFormat12(WaasDescription)
}

final case class WaasDescription(valid: Boolean,
                                 validWorkflow: Boolean,
                                 errors: List[String],
                                 name: String,
                                 inputs: Option[List[WaasDescriptionInputDefinition]],
                                 outputs: Option[List[WaasDescriptionOutputDefinition]],
                                 images: Option[List[String]],
                                 submittedDescriptorType: Option[WaasWorkflowDescriptorType],
                                 importedDescriptorTypes: Option[List[WaasWorkflowDescriptorType]],
                                 meta: Option[JsObject],
                                 parameterMeta: Option[JsObject],
                                 isRunnableWorkflow: Option[Boolean])

final case class WaasDescriptionInputDefinition(name: String,
                                                valueType: WaasDescriptionWomType,
                                                optional: Option[Boolean],
                                                default: Option[JsValue],
                                                typeDisplayName: String)

final case class WaasDescriptionOutputDefinition(name: String,
                                                 valueType: WaasDescriptionWomType,
                                                 typeDisplayName: String)

final case class WaasDescriptionWomType(typeName: String)
final case class WaasWorkflowDescriptorType(descriptorType: String, descriptorTypeVersion: String)
