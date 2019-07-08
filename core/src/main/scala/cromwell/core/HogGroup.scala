package cromwell.core

import com.typesafe.config.ConfigFactory

import scala.util.{Failure, Success}

final case class HogGroup(value: String) extends AnyVal

object HogGroup {

  type HogGroupDeciderFunction = (WorkflowOptions, WorkflowId) => HogGroup

  // NB: This is separated out from the apply so that we only have to load the config once:
  val HogGroupDeciderFunction: HogGroupDeciderFunction = {
    val config = ConfigFactory.load

    if (config.hasPath("system.hog-safety.workflow-option")) {
      val hogGroupField = config.getString("system.hog-safety.workflow-option")

      (options, workflowId) => {
        options.get(hogGroupField) match {
          case Success(hg) => HogGroup(hg)
          case Failure(_) => HogGroup(workflowId.shortString)
        }
      }
    } else {
      (_, workflowId) => HogGroup(workflowId.shortString)
    }
  }

  def decide(workflowOptions: WorkflowOptions, workflowId: WorkflowId): HogGroup = HogGroupDeciderFunction.apply(workflowOptions, workflowId)
}
