package cromiam.sam

import akka.actor.ActorSystem
import akka.event.LoggingAdapter
import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken}
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream.ActorMaterializer
import akka.util.ByteString
import com.softwaremill.sttp._
import cromiam.auth.{Collection, User}
import cromiam.sam.SamClient._
import cromiam.sam.SamResourceJsonSupport._
import cromiam.server.status.StatusCheckedSubsystem
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode

import scala.concurrent.{ExecutionContextExecutor, Future}

/*
  TODO: There exists a swagger codegen Sam client somewhere, and there also exists a Scala wrapper for it in Leo.
     I've been told that the former isn't *quite* ready for primetime and I've requested that the latter be pulled
     out into workbench-libs. We should replace this with that once the stars line up but for now it doesn't seem
     worth it. If one finds themselves making heavy changes to this file, that statement should be reevaluated.
 */
class SamClient(scheme: String,
                interface: String,
                port: Int,
                whitelistAuthModeOption: Option[GoogleAuthMode],
                log: LoggingAdapter)
               (implicit system: ActorSystem,
                ece: ExecutionContextExecutor,
                materializer: ActorMaterializer) extends StatusCheckedSubsystem {
  override val statusUri = uri"$samBaseUri/status"

  /**
    * Returns a future boolean if the user is on the whitelist, otherwise returns false.
    *
    * If the whitelist is not enabled defaults to true.
    */
  def isWhitelisted(user: User): Future[Boolean] = {
    whitelistAuthModeOption match {
      case None => Future.successful(true)
      case Some(whiteListAuthMode) =>
        // TODO: Ask Sam for the whitelist OR ask Sam if the user on the whitelist. For now just ensure we can get a token.
        // TODO: Update cromiam.server.config.WhitelistAuthMode.WhitelistScopes with the correct scopes!
        freshWhitelistAuthorization(whiteListAuthMode)
        //val request = HttpRequest(
        //  method = HttpMethods.GET,
        //  uri = samBaseUriForWhitelistUser(user),
        //  headers = List[HttpHeader](freshWhitelistAuthorization()))

        //for {
        //  response <- Http().singleRequest(request)
        //} yield ...

        Future.successful(true)
    }
  }

  def collectionsForUser(user: User): Future[List[Collection]] = {
    val request = HttpRequest(method = HttpMethods.GET, uri = samBaseCollectionUri, headers = List[HttpHeader](user.authorization))

    for {
      response <- Http().singleRequest(request)
      resources <- Unmarshal(response.entity).to[List[SamResource]]
    } yield resources.map(r => Collection(r.resourceId))
  }

  /**
    * Requests authorization for the authenticated user to perform an action from the Sam service for a collection
    * @return Successful future if the auth is accepted, a Failure otherwise.
    */
  def requestAuth(authorizationRequest: CollectionAuthorizationRequest): Future[Unit] = {
    val logString = authorizationRequest.action + " access for user " + authorizationRequest.user.userId +
      " on a request to " + authorizationRequest.action +  " for collection " + authorizationRequest.collection.name

    def validateEntityBytes(byteString: ByteString): Future[Unit] = {
      if (byteString.utf8String == "true") Future.successful(())
      else {
        log.warning("Sam denied " + logString)
        Future.failed(SamDenialException)
      }
    }

    log.info("Requesting authorization for " + logString)

    val request = HttpRequest(method = HttpMethods.GET,
      uri = samAuthorizeActionUri(authorizationRequest),
      headers = List[HttpHeader](authorizationRequest.user.authorization))

    for {
      response <- Http().singleRequest(request)
      entityBytes <- response.entity.dataBytes.runFold(ByteString(""))(_ ++ _)
      _ <- validateEntityBytes(entityBytes)
    } yield ()
  }

  /**
      - Try to create the collection
        - If 409 is the response, it already exists - check to see if user has 'add' permission
        - else we're good
      - If the result of the above is ok then we're ok.
   */
  def requestSubmission(user: User, collection: Collection): Future[Unit] = {
    log.info("Verifying user " + user.userId + " can submit a workflow to collection " + collection.name)
    val createCollection = registerCreation(user, collection)

    createCollection flatMap {
      case r if r.status == StatusCodes.Conflict => requestAuth(CollectionAuthorizationRequest(user, collection, "add"))
      case _ => Future.successful(())
    }
  }

  private def registerCreation(user: User, collection: Collection): Future[HttpResponse] = {
    val request = HttpRequest(method = HttpMethods.POST, uri = samRegisterUri(collection), headers = List[HttpHeader](user.authorization))
    Http().singleRequest(request)
  }

  private def samAuthorizeActionUri(authorizationRequest: CollectionAuthorizationRequest) = {
    akka.http.scaladsl.model.Uri(s"${samBaseUriForWorkflow(authorizationRequest.collection)}/action/${authorizationRequest.action}")
  }

  private def freshWhitelistAuthorization(whitelistAuthMode: GoogleAuthMode): Authorization = {
    Authorization(OAuth2BearerToken(whitelistAuthMode.freshAccessToken()))
  }

  private def samRegisterUri(collection: Collection) = akka.http.scaladsl.model.Uri(samBaseUriForWorkflow(collection))

  private def samBaseUriForWorkflow(collection: Collection) = s"$samBaseCollectionUri/${collection.name}"

  private lazy val samBaseUri = s"$scheme://$interface:$port"
  private lazy val samBaseResourceUri = s"$samBaseUri/api/resource"
  private lazy val samBaseCollectionUri = s"$samBaseResourceUri/workflow-collection"

}

object SamClient {
  case object SamDenialException extends Exception("Access Denied")
  val SamDenialResponse = HttpResponse(status = StatusCodes.Unauthorized, entity = SamDenialException.getMessage)
  final case class SamConnectionFailure(phase: String, f: Throwable) extends Exception(s"Unable to connect to Sam during $phase (${f.getMessage})", f)

  final case class CollectionAuthorizationRequest(user: User, collection: Collection, action: String)

}
