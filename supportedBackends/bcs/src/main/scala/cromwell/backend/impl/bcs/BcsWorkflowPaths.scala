package cromwell.backend.impl.bcs

import com.typesafe.config.Config
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.path.{PathBuilder}


case class BcsWorkflowPaths(override val workflowDescriptor: BackendWorkflowDescriptor,
                            override val config: Config,
                            override val pathBuilders: List[PathBuilder] = WorkflowPaths.DefaultPathBuilders) extends WorkflowPaths {


    override def toJobPaths(workflowPaths: WorkflowPaths, jobKey: BackendJobDescriptorKey): BcsJobPaths = {
        new BcsJobPaths(workflowPaths.asInstanceOf[BcsWorkflowPaths], jobKey)
    }

    override protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths = this.copy(workflowDescriptor = workflowDescriptor)

    private[bcs] def getWorkflowInputMounts: BcsInputMount = {
        val src = workflowRoot
        val dest = BcsJobPaths.BcsTempInputDirectory.resolve(src.pathWithoutScheme)
        BcsInputMount(src, dest, true)
    }
}
