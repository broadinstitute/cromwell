package wdl4s.wdl

trait WorkflowScoped extends Scope {
  def parentWorkflow: WdlWorkflow = ancestry.collectFirst({ case w: WdlWorkflow => w }).getOrElse(
    throw new IllegalStateException(s"Grammar constraint violation: $fullyQualifiedName should be contained in a workflow")
  )
}
