package cromwell.api.model

import spray.json.{DefaultJsonProtocol, JsObject, JsValue}
import spray.json.RootJsonFormat

object WorkflowDescriptionJsonSupport extends DefaultJsonProtocol {
  implicit val WaasWorkflowDescriptorTypeFormat: RootJsonFormat[WaasWorkflowDescriptorType] =
    jsonFormat2(WaasWorkflowDescriptorType)
  implicit val WaasDescriptionWomTypeFormat: RootJsonFormat[WaasDescriptionWomType] =
    jsonFormat1(WaasDescriptionWomType)
  implicit val WaasDescriptionInputDefinitionFormat: RootJsonFormat[WaasDescriptionInputDefinition] =
    jsonFormat5(WaasDescriptionInputDefinition)
  implicit val WaasDescriptionOutputDefinitionFormat: RootJsonFormat[WaasDescriptionOutputDefinition] =
    jsonFormat3(WaasDescriptionOutputDefinition)

  implicit val WorkflowDescriptionFormat: RootJsonFormat[WaasDescription] = jsonFormat12(WaasDescription)
}

final case class WaasDescription(valid: Boolean,
                                 validWorkflow: Boolean,
                                 errors: List[String],
                                 name: String,
                                 inputs: List[WaasDescriptionInputDefinition],
                                 outputs: List[WaasDescriptionOutputDefinition],
                                 images: List[String],
                                 submittedDescriptorType: WaasWorkflowDescriptorType,
                                 importedDescriptorTypes: List[WaasWorkflowDescriptorType],
                                 meta: JsObject,
                                 parameterMeta: JsObject,
                                 isRunnableWorkflow: Boolean
)

final case class WaasDescriptionInputDefinition(name: String,
                                                valueType: WaasDescriptionWomType,
                                                optional: Option[Boolean],
                                                default: Option[JsValue],
                                                typeDisplayName: String
)

final case class WaasDescriptionOutputDefinition(name: String,
                                                 valueType: WaasDescriptionWomType,
                                                 typeDisplayName: String
)

final case class WaasDescriptionWomType(typeName: String)
final case class WaasWorkflowDescriptorType(descriptorType: Option[String], descriptorTypeVersion: Option[String])
