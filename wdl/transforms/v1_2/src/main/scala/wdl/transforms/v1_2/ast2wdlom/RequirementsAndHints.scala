package wdl.transforms.v1_2.ast2wdlom

import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{KvPair, ObjectLiteral, StringLiteral}
import wdl.model.draft3.elements.RuntimeAttributesSectionElement
import wdl.model.v1_2.elements.{HintsSectionElement, RequirementsSectionElement}
import wom.RuntimeAttributesKeys

/**
 * Converts the WDL 1.2 `requirements` and `hints` sections into the runtime attributes understood by Cromwell's
 * backends.
 *
 * https://github.com/openwdl/wdl/blob/v1.2.1/SPEC.md#-requirements-section
 * https://github.com/openwdl/wdl/blob/v1.2.1/SPEC.md#-hints-section
 */
object RequirementsAndHints {

  // WDL 1.2 requirement name (and alias) -> Cromwell runtime attribute key
  private val requirementAttributes: Map[String, String] = Map(
    "container" -> RuntimeAttributesKeys.ContainerKey,
    "docker" -> RuntimeAttributesKeys.DockerKey,
    "cpu" -> RuntimeAttributesKeys.CpuKey,
    "memory" -> RuntimeAttributesKeys.MemoryKey,
    "gpu" -> RuntimeAttributesKeys.GpuRequiredKey,
    // TODO: No Cromwell backend can provision FPGAs. Per the spec, a task whose requirements cannot be met must fail.
    //  Suggestion: add an `fpga` validation to `StandardValidatedRuntimeAttributesBuilder` that fails the job when the
    //  value is `true`. Until then the attribute is passed through and backends ignore it with a warning.
    "fpga" -> "fpga",
    // TODO: WDL 1.2 `disks` accepts an Int (GiB), a String ("10 GiB", "/mnt/outputs 10 GiB") or an Array[String].
    //  Cromwell backends expect their own formats (e.g. "local-disk 10 SSD" on GCP Batch, "local-disk 10 SSD" or
    //  mount points on AWS). Suggestion: teach each backend's `disks` validation to also accept the WDL 1.2 formats.
    "disks" -> "disks",
    "max_retries" -> RuntimeAttributesKeys.MaxRetriesKey,
    "maxRetries" -> RuntimeAttributesKeys.MaxRetriesKey,
    "return_codes" -> RuntimeAttributesKeys.ReturnCodesKey,
    "returnCodes" -> RuntimeAttributesKeys.ReturnCodesKey
  )

  // Reserved WDL 1.2 task hint name (and alias) -> Cromwell runtime attribute key
  private val reservedHintAttributes: Map[String, String] = Map(
    "localization_optional" -> RuntimeAttributesKeys.LocalizationOptional,
    "localizationOptional" -> RuntimeAttributesKeys.LocalizationOptional,
    // TODO: Keys of `input` hints may use dotted notation for struct members (e.g. `person.cv`), but backends only look
    //  up top-level input names in `runtime.inputs` (see `StandardAsyncExecutionActor.runtimeInputsToNotLocalize`).
    //  Suggestion: resolve dotted keys to the files inside the struct input when deciding what not to localize.
    "inputs" -> RuntimeAttributesKeys.Inputs,
    "outputs" -> RuntimeAttributesKeys.Outputs,
    // TODO: `max_cpu`, `max_memory` and `short_task` are passed to the backends under their WDL 1.1 names, which no
    //  backend currently acts on. Suggestion: map `short_task: true` to `preemptible` on GCP Batch, and use
    //  `max_cpu`/`max_memory` as upper bounds when a backend picks a machine type.
    "max_cpu" -> "maxCpu",
    "maxCpu" -> "maxCpu",
    "max_memory" -> "maxMemory",
    "maxMemory" -> "maxMemory",
    "short_task" -> "shortTask",
    "shortTask" -> "shortTask"
  )

  // TODO: The reserved `disks`, `gpu` and `fpga` hints have different types than the requirements of the same name
  //  (e.g. `gpu` is a count or a specification, not a Boolean), so they cannot be passed through under the same runtime
  //  attribute key. They are ignored for now, which the spec allows. Suggestion: map an Int `gpu` hint to `gpuCount`,
  //  and a String or Map `disks` hint to additional mount points, backend by backend.
  private val ignoredHints: Set[String] = Set("disks", "gpu", "fpga")

  // Hints that may be specified inside a `hints { ... }` object (e.g. `inputs: input { foo: hints { ... } }`)
  private val nestedHintAttributes: Map[String, String] = Map(
    "localization_optional" -> RuntimeAttributesKeys.LocalizationOptional
  )

  // The deprecated `runtime` section accepts arbitrary attributes, but the `requirements` attributes are reserved and
  // have the same semantics there. Rename the WDL 1.2 spellings to the keys Cromwell understands.
  private val runtimeAttributeRenames: Map[String, String] = Map(
    "max_retries" -> RuntimeAttributesKeys.MaxRetriesKey,
    "return_codes" -> RuntimeAttributesKeys.ReturnCodesKey,
    "localization_optional" -> RuntimeAttributesKeys.LocalizationOptional
  )

  def normalizeRuntimeSection(runtime: RuntimeAttributesSectionElement): ErrorOr[RuntimeAttributesSectionElement] = {
    val renamed = runtime.runtimeAttributes.map(kv => kv.copy(key = runtimeAttributeRenames.getOrElse(kv.key, kv.key)))
    validateNoDuplicates(renamed, "runtime").map(RuntimeAttributesSectionElement.apply)
  }

  def toRuntimeAttributesSection(requirements: Option[RequirementsSectionElement],
                                 hints: Option[HintsSectionElement]
  ): ErrorOr[Option[RuntimeAttributesSectionElement]] = {
    val requirementsValidation: ErrorOr[Vector[KvPair]] =
      requirements.map(_.requirements).getOrElse(Vector.empty) traverse {
        case KvPair(key, value) if requirementAttributes.contains(key) =>
          KvPair(requirementAttributes(key), value).validNel
        case KvPair(key, _) =>
          (s"Unsupported attribute '$key' in the requirements section. Allowed attributes are: " +
            s"${requirementAttributes.keys.toList.sorted.mkString(", ")}. Other attributes belong in the hints section.").invalidNel
      }

    val hintAttributes: Vector[KvPair] = hints.map(_.hints).getOrElse(Vector.empty) collect {
      case KvPair(key, value) if reservedHintAttributes.contains(key) =>
        KvPair(reservedHintAttributes(key), normalizeNestedHints(value))
      // TODO: Hints are optional by definition: "a task execution never fails due to the inability of the execution
      //  engine to recognize or satisfy a hint". Non-reserved hints are passed through so that Cromwell-specific
      //  attributes (e.g. `preemptible`, `bootDiskSizeGb`, `zones`, `noAddress`, `queueArn`) keep working in WDL 1.2,
      //  but an invalid value will still fail the job in backend validation. Suggestion: decide whether
      //  Cromwell-specific attributes should be namespaced (e.g. `gcp: hints { preemptible: 3 }`, see "Compute
      //  Environments" in the spec), and make backend validation of hint-sourced attributes warn instead of fail.
      case KvPair(key, value) if !ignoredHints.contains(key) => KvPair(key, value)
    }

    requirementsValidation andThen { requirementAttributes =>
      validateNoDuplicates(requirementAttributes, "requirements") andThen { validRequirements =>
        // Requirements take precedence over hints of the same name
        val requirementKeys = validRequirements.map(_.key).toSet
        val allAttributes =
          (validRequirements ++ hintAttributes.filterNot(hint => requirementKeys.contains(hint.key)))
            .filterNot(isAnyContainer)

        validateContainerAliases(allAttributes) map { _ =>
          if (requirements.isEmpty && hints.isEmpty) None
          else Option(RuntimeAttributesSectionElement(allAttributes))
        }
      }
    }
  }

  // `container: "*"` means the task may run in any POSIX environment, which is what Cromwell does without a container.
  // TODO: `"*"` is only recognized as a literal. An expression evaluating to `"*"`, or an array of images containing
  //  `"*"`, is passed to the backend as an image name. Suggestion: handle `"*"` in `Containers`/`DockerValidation`.
  private def isAnyContainer(kvPair: KvPair): Boolean =
    (kvPair.key == RuntimeAttributesKeys.ContainerKey || kvPair.key == RuntimeAttributesKeys.DockerKey) &&
      kvPair.value == StringLiteral("*")

  private def validateContainerAliases(attributes: Vector[KvPair]): ErrorOr[Unit] = {
    val keys = attributes.map(_.key).toSet
    if (keys.contains(RuntimeAttributesKeys.ContainerKey) && keys.contains(RuntimeAttributesKeys.DockerKey))
      "The requirements section cannot specify both 'container' and its alias 'docker'".invalidNel
    else ().validNel
  }

  private def validateNoDuplicates(attributes: Vector[KvPair], sectionName: String): ErrorOr[Vector[KvPair]] = {
    val duplicates = attributes.groupBy(_.key).collect { case (key, kvs) if kvs.size > 1 => key }
    if (duplicates.isEmpty) attributes.validNel
    else
      s"Duplicate attributes in the $sectionName section (including aliases): ${duplicates.toList.sorted.mkString(", ")}".invalidNel
  }

  // Rename reserved hints nested inside `input`/`output`/`hints` objects, e.g.
  //   inputs: input { foo: hints { localization_optional: true } }
  // becomes the WDL 1.1-style runtime attribute Cromwell understands:
  //   inputs: object { foo: object { localizationOptional: true } }
  private def normalizeNestedHints(value: ExpressionElement): ExpressionElement = value match {
    case ObjectLiteral(elements) =>
      ObjectLiteral(elements map { case (key, nested) =>
        nestedHintAttributes.getOrElse(key, key) -> normalizeNestedHints(nested)
      })
    case other => other
  }
}
