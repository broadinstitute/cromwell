package cromwell.engine.workflow

import akka.actor.{Actor, ActorRef, Props}
import akka.pattern.{ask, pipe}
import com.typesafe.scalalogging.LazyLogging
import cromwell.backend.BackendValidationActor.{FailedValidationResult, ValidationResult, Validate}
import cromwell.backend.DefaultBackendFactory
import cromwell.backend.config.{BackendConfiguration, BackendConfigurationEntry}
import cromwell.engine.CromwellActor
import cromwell.webservice.APIResponse
import cromwell.webservice.PerRequest.RequestComplete
import spray.http.StatusCodes
import spray.json._
import wdl4s._
import spray.httpx.SprayJsonSupport._
import cromwell.webservice.WorkflowJsonSupport._

import scala.concurrent.duration.Duration
import scala.concurrent.{Promise, Await, Future}
import scala.language.postfixOps
import scala.util.{Try, Failure, Success}


object ValidateActor {
  private val tag = "ValidateActor"

  def props(wdlSource: WdlSource, wdlJson: Option[WdlJson], workflowOptions: Option[WdlJson]): Props = {
    Props(new ValidateActor(wdlSource, wdlJson, workflowOptions))
  }

  sealed trait ValidateActorMessage

  //TODO: [gaurav] Check what and how is this being used
  implicit class EnhancedCall(val call: Call) extends AnyVal {
    def toRuntimeAttributes = call.task.runtimeAttributes
  }

  case object ValidateWorkflow extends ValidateActorMessage

  case object GetExecutableBackends extends ValidateActorMessage

}

class ValidateActor(wdlSource: WdlSource, workflowInputs: Option[WdlJson], workflowOptions: Option[String])
  extends Actor with CromwellActor with LazyLogging {

  import ValidateActor._
  import context.dispatcher

  override def receive = {
    case ValidateWorkflow =>
      validateWorkflow(sender())
    // NOTE: self shuts down when the parent PerRequest shuts down
    case GetExecutableBackends =>
      val requster = sender()
      val futureValidatedBackends: Future[Seq[BackendConfigurationEntry]] = for {
        namespaceWithWorkflow <- Future(NamespaceWithWorkflow.load(wdlSource))
        inputs <- Future(workflowInputs.get.parseJson).map(_.asJsObject.fields)
        coercedInputs <- Future.fromTry(namespaceWithWorkflow.coerceRawInputs(inputs))
        validatedBackends <- executableBackendsForWf(namespaceWithWorkflow, Option(coercedInputs), workflowOptions)
      } yield (validatedBackends)

      val x = Await.result(futureValidatedBackends, Duration.Inf)

      futureValidatedBackends pipeTo requster
  }

  //TODO: This should be removed / deprecated?
  private def validateWorkflow(sentBy: ActorRef): Unit = {
    logger.info(s"$tag for $sentBy")
    val futureValidation: Future[Unit] = for {
      namespaceWithWorkflow <- Future(NamespaceWithWorkflow.load(wdlSource))
      inputs <- Future(workflowInputs.get.parseJson).map(_.asJsObject.fields)
      coercedInputs <- Future.fromTry(namespaceWithWorkflow.coerceRawInputs(inputs))
      runtime = namespaceWithWorkflow.workflow.calls foreach {
        _.toRuntimeAttributes
      }
    } yield () // Validate that the future run and return `Success[Unit]` aka (), or `Failure[Exception]`

    futureValidation onComplete {
      case Success(_) =>
        logger.info(s"$tag success $sentBy")
        sentBy ! RequestComplete(
          StatusCodes.OK,
          APIResponse.success("Validation succeeded."))

      case Failure(ex) =>
        val messageOrBlank = Option(ex.getMessage).mkString
        logger.info(s"$tag error $sentBy: $messageOrBlank")
        sentBy ! RequestComplete(
          StatusCodes.BadRequest,
          APIResponse.fail(ex))
    }
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

    if (!willExecute.isEmpty) logger.info(s"Backends with no validation requirements: $willExecute")

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
        case (k, v: FailedValidationResult) => logger.info(s"Backend [${k.name} failed validation. Errors: ${v.errors}")
        case oops => logger.error(s"Oops! This is embarrasingly unexpected! [$oops]")
      }
    }

    // This will combine the Successfully validated backends with the list `willExecute`
    val aggregatedValidatedBackends = for {
      backendEntry <- validatedBackends.map(_.keysIterator.toList)
      aggregatedBackendSeq = backendEntry flatMap (backend => willExecute.:+(backend)) toSeq
    } yield (aggregatedBackendSeq)

    //This `promise` is here because the above for comprehension may not combine the two lists
    // i.e. the successful `mayExecute` and `willExecute` if the former is empty. As such, we need to
    //modify the value of that future to just return the `willExecute` list
    val promise = Promise[Seq[BackendConfigurationEntry]]()

    aggregatedValidatedBackends onComplete {
      case Success(listOfBackends) =>
        if (listOfBackends.isEmpty) promise.complete(Try(willExecute)) else promise.complete(Try(listOfBackends))
      case Failure(reason) =>
        logger.error(s"Something went wrong while getting validated backends! Msg: ${reason.getMessage}", reason)
        promise.complete(Try(Seq()))
    }
    promise.future
  }
}
