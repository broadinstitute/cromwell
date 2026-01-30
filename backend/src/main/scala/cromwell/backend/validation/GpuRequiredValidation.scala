package cromwell.backend.validation

import wom.RuntimeAttributesKeys
import wom.values.WomBoolean

/**
 * This runtime attribute indicates whether GPU resources are required for the job. If set to true, the backend
 * must ensure that the execution environment has access to GPU resources; if not, the task should fail.
 * Supported starting in WDL 1.1
 * https://github.com/openwdl/wdl/blob/wdl-1.1/SPEC.md#gpu
 */

object GpuRequiredValidation extends BooleanRuntimeAttributesValidation(RuntimeAttributesKeys.GpuRequiredKey) {
  val DefaultValue: WomBoolean = WomBoolean(false)
}
