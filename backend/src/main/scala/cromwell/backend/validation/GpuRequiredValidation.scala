package cromwell.backend.validation

import wom.RuntimeAttributesKeys

/**
 * This runtime attribute indicates whether GPU resources are required for the job. If set to true, the backend
 * must ensure that the execution environment has access to GPU resources; if not, the task should fail.
 * Supported starting in WDL 1.1
 * https://github.com/openwdl/wdl/blob/wdl-1.1/SPEC.md#gpu
 */

object GpuRequiredValidation {
  lazy val instance: RuntimeAttributesValidation[Boolean] = new GpuRequiredValidation
  lazy val optional: OptionalRuntimeAttributesValidation[Boolean] = instance.optional
}

class GpuRequiredValidation extends BooleanRuntimeAttributesValidation(RuntimeAttributesKeys.GpuRequiredKey)
