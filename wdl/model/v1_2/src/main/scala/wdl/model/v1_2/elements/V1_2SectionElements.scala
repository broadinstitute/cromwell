package wdl.model.v1_2.elements

import wdl.model.draft3.elements.ExpressionElement.KvPair
import wdl.model.draft3.elements.{TaskSectionElement, WorkflowBodyElement}
import wom.callable.MetaKvPair

/** A task `requirements` section (WDL 1.2). Only the attributes defined by the specification are allowed. */
final case class RequirementsSectionElement(requirements: Vector[KvPair]) extends TaskSectionElement

/**
 * A task `hints` section (WDL 1.2). Any attribute is allowed. The hints-scoped `hints`, `input` and `output` literals
 * are represented as object literals.
 */
final case class HintsSectionElement(hints: Vector[KvPair]) extends TaskSectionElement

/** A workflow `hints` section (WDL 1.2). Workflow hints must have literal values. */
final case class WorkflowHintsSectionElement(hints: Vector[MetaKvPair]) extends WorkflowBodyElement
