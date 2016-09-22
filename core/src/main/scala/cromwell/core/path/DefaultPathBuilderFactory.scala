package cromwell.core.path

import cromwell.core.WorkflowOptions

case object DefaultPathBuilderFactory extends PathBuilderFactory {
  override def withOptions(options: WorkflowOptions) = DefaultPathBuilder
}
