package cwl.requirement

import cwl._
import shapeless.Poly1
import wom.RuntimeAttributesKeys._
import wom.expression.{ValueAsAnExpression, WomExpression}
import wom.values.WomString

object RequirementToAttributeMap extends Poly1 {
  implicit def fromJs: Case.Aux[InlineJavascriptRequirement, Set[String] => Map[String, WomExpression]] = at[InlineJavascriptRequirement] {
    _ => _ => Map.empty
  }

  implicit def fromSchemaDef: Case.Aux[SchemaDefRequirement, Set[String] => Map[String, WomExpression]] = at[SchemaDefRequirement] {
    _ => _ => Map.empty
  }

  implicit def fromDocker: Case.Aux[DockerRequirement, Set[String] => Map[String, WomExpression]] = at[DockerRequirement] {
    docker => _ => docker.dockerPull.orElse(docker.dockerImageId).map({ pull =>
      DockerKey -> ValueAsAnExpression(WomString(pull))
    }).toMap
  }

  implicit def fromSoftware: Case.Aux[SoftwareRequirement, Set[String] => Map[String, WomExpression]] = at[SoftwareRequirement] {
    _ => _ => Map.empty
  }

  implicit def fromInitialWorkDir: Case.Aux[InitialWorkDirRequirement, Set[String] => Map[String, WomExpression]] = at[InitialWorkDirRequirement] {
    _ => _ => Map.empty
  }

  implicit def fromEnvVar: Case.Aux[EnvVarRequirement, Set[String] => Map[String, WomExpression]] = at[EnvVarRequirement] {
    _ => _ => Map.empty
  }

  implicit def fromShellCommand: Case.Aux[ShellCommandRequirement, Set[String] => Map[String, WomExpression]] = at[ShellCommandRequirement] {
    _ => _ => Map.empty
  }

  implicit def fromResource: Case.Aux[ResourceRequirement, Set[String] => Map[String, WomExpression]] = at[ResourceRequirement] {
    resource => inputNames =>
      def toExpression(resourceRequirement: ResourceRequirementType) = resourceRequirement.fold(ResourceRequirementToWomExpression).apply(inputNames)

      List(
        // Map cpuMin to both cpuMin and cpu keys
        resource.effectiveCoreMin.toList.map(toExpression).flatMap(min => List(CpuMinKey -> min, CpuKey -> min)),
        resource.effectiveCoreMax.toList.map(toExpression).map(CpuMaxKey -> _),
        // Map ramMin to both memoryMin and memory keys
        resource.effectiveRamMin.toList.map(toExpression).flatMap(min => List(MemoryMinKey -> min, MemoryKey -> min)),
        resource.effectiveRamMax.toList.map(toExpression).map(MemoryMaxKey -> _),
        resource.effectiveTmpdirMin.toList.map(toExpression).map(TmpDirMinKey -> _),
        resource.effectiveTmpdirMax.toList.map(toExpression).map(TmpDirMaxKey -> _),
        resource.effectiveOutdirMin.toList.map(toExpression).map(OutDirMinKey -> _),
        resource.effectiveOutdirMax.toList.map(toExpression).map(OutDirMaxKey -> _)
      ).flatten.toMap
  }

  implicit def fromSubWorkflow: Case.Aux[SubworkflowFeatureRequirement, Set[String] => Map[String, WomExpression]] = at[SubworkflowFeatureRequirement] {
    _ => _ => Map.empty
  }

  implicit def fromScatter: Case.Aux[ScatterFeatureRequirement, Set[String] => Map[String, WomExpression]] = at[ScatterFeatureRequirement] {
    _ => _ => Map.empty
  }

  implicit def fromMultipleInput: Case.Aux[MultipleInputFeatureRequirement, Set[String] => Map[String, WomExpression]] = at[MultipleInputFeatureRequirement] {
    _ => _ => Map.empty
  }

  implicit def fromStepInput: Case.Aux[StepInputExpressionRequirement, Set[String] => Map[String, WomExpression]] = at[StepInputExpressionRequirement] {
    _ => _ => Map.empty
  }
}
