package cromwell.backend.impl.sfs.config

import cromwell.backend.validation.{Containers, RuntimeAttributesValidation, ValidatedRuntimeAttributes}
import wdl.draft2.model.Declaration
import wom.values.{WomString, WomValue}

class ContainerDeclarationValidation(declaration: Declaration, instanceValidation: RuntimeAttributesValidation[_])
    extends DeclarationValidation(declaration, instanceValidation, usedInCallCachingOverride = None) {

  override def extractWdlValueOption(validatedRuntimeAttributes: ValidatedRuntimeAttributes): Option[WomValue] =
    RuntimeAttributesValidation.extractOption[Containers](instanceValidation.key, validatedRuntimeAttributes) map {
      containers => WomString(containers.values.head)
    }
}
