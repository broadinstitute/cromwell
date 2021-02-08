package cromwell.backend.impl.sfs.config

import cromwell.backend.validation.{RuntimeAttributesValidation, ValidatedRuntimeAttributes}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import wdl.draft2.model.Declaration
import wom.values.WomValue

class CpuDeclarationValidation(declaration: Declaration, instanceValidation: RuntimeAttributesValidation[Int Refined Positive])
  extends DeclarationValidation(declaration, instanceValidation, usedInCallCachingOverride = Option(false)) {

  override def extractWdlValueOption(validatedRuntimeAttributes: ValidatedRuntimeAttributes): Option[WomValue] = {
    RuntimeAttributesValidation.extractOption(instanceValidation, validatedRuntimeAttributes) map { refined =>
      declaration.womType.coerceRawValue(refined.value).get
    }
  }
}
