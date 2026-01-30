package cromwell.backend.validation

import wom.RuntimeAttributesKeys

/**
  * Different WDL versions support different names for the runtime attribute that specifies the container image to use.
  * WDL 1.0 supports only `docker`, WDL 1.1 and later support `container` (preferred) and `docker` (deprecated).
  */
object DockerValidation {
  lazy val instance: OptionalRuntimeAttributesValidation[Containers] = new DockerValidation
}

class DockerValidation extends ContainersValidation {
  override val key: String = RuntimeAttributesKeys.DockerKey
}

object ContainerValidation {
  lazy val instance: OptionalRuntimeAttributesValidation[Containers] = new ContainerValidation
}

class ContainerValidation extends ContainersValidation {
  override val key: String = RuntimeAttributesKeys.ContainerKey
}
