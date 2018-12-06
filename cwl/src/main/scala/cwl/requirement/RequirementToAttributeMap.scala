package cwl.requirement

import cwl._
import shapeless.Poly1
import wom.RuntimeAttributesKeys._
import wom.expression.{ValueAsAnExpression, WomExpression}
import wom.values.{WomLong, WomString}

object RequirementToAttributeMap extends Poly1 {
  type ResourcesToExpressionMap = (Set[String], ExpressionLib) => Map[String, WomExpression]
  implicit def fromJs: Case.Aux[InlineJavascriptRequirement, ResourcesToExpressionMap] = at[InlineJavascriptRequirement] {
    _ => (_,_) => Map.empty
  }

  implicit def fromSchemaDef: Case.Aux[SchemaDefRequirement, ResourcesToExpressionMap] = at[SchemaDefRequirement] {
    _ => (_,_) => Map.empty
  }

  implicit def fromDocker: Case.Aux[DockerRequirement, ResourcesToExpressionMap] = at[DockerRequirement] {
    docker => (_,_) => docker.dockerPull.orElse(docker.dockerImageId).map({ pull =>
      DockerKey -> ValueAsAnExpression(WomString(pull))
    }).toMap
  }

  implicit def fromSoftware: Case.Aux[SoftwareRequirement, ResourcesToExpressionMap] = at[SoftwareRequirement] {
    _ => (_,_) => Map.empty
  }

  implicit def fromInitialWorkDir: Case.Aux[InitialWorkDirRequirement, ResourcesToExpressionMap] = at[InitialWorkDirRequirement] {
    _ => (_,_) => Map.empty
  }

  implicit def fromEnvVar: Case.Aux[EnvVarRequirement, ResourcesToExpressionMap] = at[EnvVarRequirement] {
    _ => (_,_) => Map.empty
  }

  implicit def fromShellCommand: Case.Aux[ShellCommandRequirement, ResourcesToExpressionMap] = at[ShellCommandRequirement] {
    _ => (_,_) => Map.empty
  }

  implicit def fromResource: Case.Aux[ResourceRequirement, ResourcesToExpressionMap] = at[ResourceRequirement] {
    resource => (inputNames, expressionLib) =>
      def toExpression(resourceRequirement: ResourceRequirementType) =
        resourceRequirement.fold(ResourceRequirementToWomExpression).apply(inputNames, expressionLib)

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

  implicit def fromInputResourceRequirement: Case.Aux[DnaNexusInputResourceRequirement, ResourcesToExpressionMap] = at[DnaNexusInputResourceRequirement] {
    case DnaNexusInputResourceRequirement(_, indirMin) => (_, _) => indirMin.map(value => ValueAsAnExpression(WomLong(value))).map(DnaNexusInputDirMinKey -> _).toMap
  }

  implicit def fromSubWorkflow: Case.Aux[SubworkflowFeatureRequirement, ResourcesToExpressionMap] = at[SubworkflowFeatureRequirement] {
    _ => (_,_) => Map.empty
  }

  implicit def fromScatter: Case.Aux[ScatterFeatureRequirement, ResourcesToExpressionMap] = at[ScatterFeatureRequirement] {
    _ => (_,_) => Map.empty
  }

  implicit def fromMultipleInput: Case.Aux[MultipleInputFeatureRequirement, ResourcesToExpressionMap] = at[MultipleInputFeatureRequirement] {
    _ => (_,_) => Map.empty
  }

  implicit def fromStepInput: Case.Aux[StepInputExpressionRequirement, ResourcesToExpressionMap] = at[StepInputExpressionRequirement] {
    _ => (_,_) => Map.empty
  }
}
