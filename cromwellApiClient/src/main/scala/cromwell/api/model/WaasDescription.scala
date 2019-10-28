package cromwell.api.model

import spray.json.{DefaultJsonProtocol, JsValue}

object WorkflowDescriptionJsonSupport extends DefaultJsonProtocol {
  implicit val WaasWorkflowDescriptorTypeFormat = jsonFormat2(WaasWorkflowDescriptorType)
  implicit val WaasDescriptionWomTypeFormat = jsonFormat1(WaasDescriptionWomType)
  implicit val WaasDescriptionParameterFormat = jsonFormat5(WaasDescriptionParameter)

  implicit val WorkflowDescriptionFormat = jsonFormat8(WaasDescription)
}

final case class WaasDescription(valid: Boolean,
                                 validWorkflow: Boolean,
                                 errors: List[String],
                                 name: String,
                                 inputs: List[WaasDescriptionParameter],
                                 outputs: List[WaasDescriptionParameter],
                                 submittedDescriptorType: WaasWorkflowDescriptorType,
                                 isRunnableWorkflow: Boolean)

final case class WaasDescriptionParameter(name: String,
                                          valueType: WaasDescriptionWomType,
                                          optional: Boolean,
                                          default: Option[JsValue],
                                          typeDisplayName: String)

final case class WaasDescriptionWomType(typeName: String)
final case class WaasWorkflowDescriptorType(descriptorType: String, descriptorTypeVersion: String)
