package cromwell

package object backend {
  /** Represents the jobKeys executed by a (potentially sub-) workflow at a given point in time */
  type JobExecutionMap = Map[BackendWorkflowDescriptor, List[BackendJobDescriptorKey]]
}
