package cromwell.backend.google.pipelines.batch

import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.api.Refined
import com.typesafe.config.Config

import cromwell.backend.validation.ValidatedRuntimeAttributes
//import cromwell.backend.google.pipelines.common.PipelinesApiConfiguration
//import cromwell.backend.google.pipelines.common.PipelinesApiRuntimeAttributes.cpuValidation
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.validation.{CpuValidation, OptionalRuntimeAttributesValidation, RuntimeAttributesValidation}
import wom.RuntimeAttributesKeys.CpuKey
//import wom.RuntimeAttributesKeys.{CpuKey, CpuMaxKey, CpuMinKey}
import wom.values.{WomInteger, WomValue}


case class GcpBatchRuntimeAttributes(cpu: Int Refined Positive)
object CpuValidation {
  lazy val instance: RuntimeAttributesValidation[Int Refined Positive] = new CpuValidation(CpuKey)
  lazy val optional: OptionalRuntimeAttributesValidation[Int Refined Positive] = instance.optional
  //lazy val instanceMin: RuntimeAttributesValidation[Int Refined Positive] = new CpuValidation(CpuMinKey)
  //lazy val optionalMin: OptionalRuntimeAttributesValidation[Int Refined Positive] = instanceMin.optional
  //lazy val instanceMax: RuntimeAttributesValidation[Int Refined Positive] = new CpuValidation(CpuMaxKey)
  //lazy val optionalMax: OptionalRuntimeAttributesValidation[Int Refined Positive] = instanceMax.optional

  lazy val defaultMin: WomValue = WomInteger(1)
  def configDefaultWomValue(config: Option[Config]): Option[WomValue] = instance.configDefaultWomValue(config)
}
object GcpBatchRuntimeAttributes {
  private def cpuValidation(runtimeConfig: Option[Config]): RuntimeAttributesValidation[Int Refined Positive] = CpuValidation
    .instance
    .withDefault(CpuValidation
      .configDefaultWomValue(runtimeConfig) getOrElse CpuValidation
      .defaultMin)

  def runtimeAttributesBuilder(batchConfiguration: GcpBatchConfiguration): StandardValidatedRuntimeAttributesBuilder = {
    val runtimeConfig = batchConfiguration.runtimeConfig
    StandardValidatedRuntimeAttributesBuilder.default(runtimeConfig).withValidation(
      cpuValidation(runtimeConfig),
    )
  }

  def apply(validatedRuntimeAttributes: ValidatedRuntimeAttributes, runtimeAttrsConfig: Option[Config]): GcpBatchRuntimeAttributes = {
    val cpu: Int Refined Positive = RuntimeAttributesValidation.extract(cpuValidation(runtimeAttrsConfig), validatedRuntimeAttributes)

    new GcpBatchRuntimeAttributes(
      cpu
    )
  }


}
