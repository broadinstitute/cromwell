package broad.cwl.model

import shapeless.{:+:, CNil}

case class Workflow(
  cwlVersion: String,
  `class`: String,
  inputs: Map[InputParameter#Id, InputParameter] :+: Map[InputParameter#Id, InputParameter#`type`] :+: Array[InputParameter] :+: CNil,
  outputs: Map[WorkflowOutputParameter#Id, WorkflowOutputParameter] :+: Map[WorkflowOutputParameter#Id, WorkflowOutputParameter#`type`] :+: Array[WorkflowOutputParameter] :+: CNil,
  steps: Map[String, WorkflowStep] :+: Array[WorkflowStep] :+: CNil
  )

