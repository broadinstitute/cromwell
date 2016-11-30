package wdl4s

trait WorkflowScoped extends Scope {
  def parentWorkflow: Workflow = ancestry.collectFirst({ case w: Workflow => w }).getOrElse(
    throw new IllegalStateException(s"Grammar constraint violation: $fullyQualifiedName should be contained in a workflow")
  )
}
