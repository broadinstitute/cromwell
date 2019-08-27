package cromwell.backend.impl.bcs

import com.typesafe.config.Config
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.path.{Path, PathBuilder}

object BcsWorkflowPaths {
    val WorkFlowTagKey = "bcs_workflow_tag"
}

case class BcsWorkflowPaths(override val workflowDescriptor: BackendWorkflowDescriptor,
                            override val config: Config,
                            override val pathBuilders: List[PathBuilder] = WorkflowPaths.DefaultPathBuilders) extends WorkflowPaths {

    import BcsWorkflowPaths._
    override def toJobPaths(workflowPaths: WorkflowPaths, jobKey: BackendJobDescriptorKey): BcsJobPaths = {
        new BcsJobPaths(workflowPaths.asInstanceOf[BcsWorkflowPaths], jobKey)
    }

    override protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths = this.copy(workflowDescriptor = workflowDescriptor)

    override protected def workflowPathBuilder(root: Path): Path = {
        workflowDescriptor.breadCrumbs.foldLeft(root)((acc, breadCrumb) => {
            breadCrumb.toPath(acc)
        }).resolve(workflowDescriptor.callable.name).resolve(tag).resolve(workflowDescriptor.id.toString + "/")
    }

    var tag: String = {
        workflowDescriptor.workflowOptions.get(WorkFlowTagKey).getOrElse("")
    }

    private[bcs] def getWorkflowInputMounts: BcsInputMount = {
        val src = workflowRoot
        val dest = BcsJobPaths.BcsTempInputDirectory.resolve(src.pathWithoutScheme)
        BcsInputMount(Left(src), Left(dest), true)
    }
}
