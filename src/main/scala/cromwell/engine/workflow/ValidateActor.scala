package cromwell.engine.workflow

import akka.actor.{Actor, ActorRef, Props}
import akka.pattern.{ask, pipe}
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.BackendValidationActor.{FailedValidationResult, ValidationResult, Validate}
import cromwell.backend.DefaultBackendFactory
import cromwell.backend.config.{BackendConfiguration, BackendConfigurationEntry}
import cromwell.engine.{WorkflowSourceFiles, CromwellActor}
import cromwell.webservice.APIResponse
import cromwell.webservice.PerRequest.RequestComplete
import spray.http.StatusCodes
import spray.json._
import wdl4s._
import spray.httpx.SprayJsonSupport._
import cromwell.webservice.WorkflowJsonSupport._

import scala.concurrent.{Promise, Future}
import scala.language.postfixOps
import scala.util.{Try, Failure, Success}

object ValidateActor {
  private val tag = "ValidateActor"

  def props(workflowSources: WorkflowSourceFiles): Props = Props(new ValidateActor(workflowSources.wdlSource, Option(workflowSources.inputsJson), Option(workflowSources.workflowOptionsJson)))

  sealed trait ValidateActorMessage
  //TODO: `WorkflowSourceFiles` should not exist here in this file! It's only here because some of the things were breaking in WF Desc if we don't have access to this
  case class ValidationSuccess(namespaceWithWorkflow: NamespaceWithWorkflow,
                               coercedInputs: Option[WorkflowCoercedInputs],
                               workflowOptions: Option[WdlJson],
                               backends: Seq[BackendConfigurationEntry],
                               workflowSources: Option[WorkflowSourceFiles] = None) extends ValidateActorMessage
  case class ValidationFailure(reason: Throwable) extends ValidateActorMessage
  case object ValidateWorkflow extends ValidateActorMessage
  case object RequestValidation extends ValidateActorMessage
}

class ValidateActor(wdlSource: WdlSource, workflowInputs: Option[WdlJson], workflowOptions: Option[String])
  extends Actor with CromwellActor with LazyLogging {

  import ValidateActor._
  import context.dispatcher

  override def receive = {
    //TODO: <b> Deprecate this. <b> It's currently being used only for the endpoint. IMO, the responsibility for returning a status code
    // should be a part of the client. This class should only return a validation result, and the callee should convert that to a
    // status code or anything that it wants
    case ValidateWorkflow =>
      val requester = sender()

      val futureValidatedBackends: Future[Seq[BackendConfigurationEntry]] = for {
        namespaceWithWorkflow <- Future(NamespaceWithWorkflow.load(wdlSource))
        inputs <- Future(workflowInputs.get.parseJson).map(_.asJsObject.fields)
        coercedInputs <- Future.fromTry(namespaceWithWorkflow.coerceRawInputs(inputs))
        validatedBackends <- executableBackendsForWf(namespaceWithWorkflow, Option(coercedInputs), workflowOptions)
      } yield validatedBackends

      futureValidatedBackends onComplete {
        case Success(validatedBacends) if validatedBacends.nonEmpty =>
          logger.info(s"$tag success $requester")
          requester ! RequestComplete(
            StatusCodes.OK,
            APIResponse.success("Validation succeeded."))

        case Success(_) =>
          val msg = "No backend available that can accept this workflow!"
          logger.info(msg)
          requester ! RequestComplete(
            StatusCodes.BadRequest,
            APIResponse.fail(new IllegalStateException(msg)))

        case Failure(ex) =>
          val messageOrBlank = Option(ex.getMessage).mkString
          logger.info(s"$tag error $requester: $messageOrBlank")
          requester ! RequestComplete(
            StatusCodes.BadRequest,
            APIResponse.fail(ex))
      }

      //This message initiates the Validation of the workflow. It will:
      // 1.) Try to build a NamespaceWithWorkflow object
      // 2.) Try to coerce the inputs
      // 3.) Find Backends that can accept this worklfow (fails this hurdle if there is no backend found)
    case RequestValidation =>
      val requster = sender()
      val futureValidationOutcome: Future[ValidationSuccess] = for {
        namespaceWithWorkflow <- Future(NamespaceWithWorkflow.load(wdlSource))
        coercedInputs <- buildCoercedInputs(namespaceWithWorkflow, workflowInputs)
        validatedBackends <- executableBackendsForWf(namespaceWithWorkflow, Option(coercedInputs), workflowOptions)
      } yield ValidationSuccess(namespaceWithWorkflow,
        Option(coercedInputs),
        workflowOptions,
        validatedBackends,
        Option(WorkflowSourceFiles(wdlSource, workflowInputs.getOrElse("{}"), workflowOptions.getOrElse("{}"))))

      futureValidationOutcome onComplete {
        case Success(validationOutcome: ValidationSuccess) =>
          requster ! validationOutcome
        case Failure(reason) =>
          requster ! ValidationFailure(reason)
      }
  }

  private def buildCoercedInputs(namespaceWithWorkflow: NamespaceWithWorkflow, rawInputs: Option[WdlJson]): Future[WorkflowCoercedInputs] = {
    val inputs = workflowInputs.get.parseJson.asJsObject.fields
    Future.fromTry(namespaceWithWorkflow.coerceRawInputs(inputs))
  }

  /**
    * This method will validate the workflow with all the available backends. The WF fails
    * validation if this method returns a list of size 0, which meant that no backend was willing to accept
    * this workflow.
    * Later on, this method can be changed to a 'per-call-basis' acceptance test by the backends,
    * i.e. it may return a Map of [Call -> ExecutableBackends]
    *
    * @param namespace     workflow with the namespaces already resolved
    * @param wfInputs      Already coerced inputs
    * @param wfOptionsJson workflow options
    * @return List of BackendConfigurationEntry objects, representing which backends can accept this WF
    */
  private def executableBackendsForWf(namespace: NamespaceWithWorkflow,
                                      wfInputs: Option[WorkflowCoercedInputs] = None,
                                      wfOptionsJson: Option[String] = None): Future[Seq[BackendConfigurationEntry]] = {
    val listOfPluggedBackends: List[BackendConfigurationEntry] = BackendConfiguration.apply().getAllBackendConfigurations()

    logger.info(s"Plugged in Backends: $listOfPluggedBackends")

    //If the validation class is not specified, it is assumed that backend will accept any WF
    // The `willExecute` list contains backends that don't have a validation class, and will accept this workflow
    val (mayExecute, willExecute) = listOfPluggedBackends.partition(_.validationClass.isDefined)

    if (willExecute.nonEmpty) logger.info(s"Backends with no validation requirements: $willExecute")

    val backendEntryValidatorActorPair: Map[BackendConfigurationEntry, ActorRef] = mayExecute map { entry =>
      entry -> DefaultBackendFactory.getBackendActorFor(entry.validationClass.get, context.system, namespace, wfInputs, wfOptionsJson)
    } toMap

    // The implicit timeout required for Ask below comes from `CromwellActor`
    val mapOfEntryToFutureValidationRes: Map[BackendConfigurationEntry, Future[ValidationResult]] = backendEntryValidatorActorPair mapValues { validateActor =>
      (validateActor ? Validate).mapTo[ValidationResult]
    }

    val futureMapOfEntryToValidationRes: Future[Map[BackendConfigurationEntry, ValidationResult]] = Future.sequence {
      mapOfEntryToFutureValidationRes map { case (entry, futureValidationResult) =>
        futureValidationResult.map((entry, _))
      }
    }.map(_.toMap)

    //Just logging
    futureMapOfEntryToValidationRes foreach {
      _ foreach { case (k, v) => logger.info(s"For Backend [${k.name}], Result: [$v]") }
    }

    //Map which does not have any backends which failed validation
    val validatedBackends: Future[Map[BackendConfigurationEntry, ValidationResult]] = futureMapOfEntryToValidationRes map {
      _.filter(!_._2.isInstanceOf[FailedValidationResult])
    }

    //Just logging
    futureMapOfEntryToValidationRes map {
      _.filter(_._2.isInstanceOf[FailedValidationResult])
    } foreach {
      _.foreach {
        case (k, v: FailedValidationResult) => logger.error(s"Backend [${k.name} failed validation. Errors: ${v.errors}")
        case oops => logger.error(s"Oops! This is embarrasingly unexpected! [$oops]")
      }
    }

    // This will combine the Successfully validated backends with the list `willExecute`
    val aggregatedValidatedBackends = for {
      backendEntry <- validatedBackends.map(_.keysIterator.toList)
      aggregatedBackendSeq = backendEntry flatMap (backend => willExecute.:+(backend)) toSeq
    } yield aggregatedBackendSeq

    //This `promise` is here because the above for comprehension may not combine the two lists
    // i.e. the successful `mayExecute` and `willExecute` if the former is empty. As such, we need to
    //modify the value of that future to just return the `willExecute` list
    val promise = Promise[Seq[BackendConfigurationEntry]]()

    aggregatedValidatedBackends onComplete {
      case Success(seqOfBackends) =>
        if (seqOfBackends.nonEmpty)
          promise.complete(Try(seqOfBackends))
        else if (willExecute.nonEmpty)
          promise.complete(Try(willExecute))
        else throw new IllegalStateException("Failed to validate this workflow against any of the available backends!")
      case Failure(reason) =>
        val msg = s"Something went wrong while getting validated backends! Msg: ${reason.getMessage}"
        logger.error(msg, reason)
        throw new IllegalStateException(msg, reason)
    }
    promise.future
  }
}
