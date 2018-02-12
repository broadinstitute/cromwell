package cromiam.webservice

import akka.http.scaladsl.server._
import akka.http.scaladsl.server.Directives._
import akka.stream.ActorMaterializer
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient
import akka.http.scaladsl.model.{HttpEntity, _}
import cromiam.auth.{Collection, User}
import akka.event.LoggingAdapter
import Collection.{CollectionLabelName, LabelContainsCollectionException}
import QuerySupport._

import scala.concurrent.ExecutionContextExecutor
import scala.util.{Failure, Success, Try}

trait QuerySupport extends RequestSupport {
  val cromwellClient: CromwellClient
  val samClient: SamClient

  val log: LoggingAdapter

  implicit def executor: ExecutionContextExecutor
  implicit val materializer: ActorMaterializer

  def queryGetRoute: Route = path("api" / "workflows" / Segment / "query") { _ =>
    get {
      preprocessQuery("GET") { (user, collections, request) =>
        processLabelsForGetQuery(user, collections) { uri =>
          val requestToForward = HttpRequest(HttpMethods.GET, uri, request.headers)
          complete {
            cromwellClient.forwardToCromwell(requestToForward)
          }
        }
      }
    }
  }

  def queryPostRoute: Route = path("api" / "workflows" / Segment / "query") { _ =>
    post {
      preprocessQuery("GET") { (user, collections, request) =>
        processLabelsForPostQuery(user, collections) { entity =>
          complete { cromwellClient.forwardToCromwell(request.withEntity(entity)) }
        }
      }
    }
  }

  /**
    * Handles shared stuff between the POST & GET Query requests. Pulls out the user information, logs the request,
    * retrieves the collections for the user, grabs the underlying HttpRequest and forwards it on to the specific
    * directive
    */
  private def preprocessQuery(method: String): Directive[(User, List[Collection], HttpRequest)] = {
    extractUser flatMap  { user =>
      log.info("Received query " + method + " request for user " + user.userId)

      onComplete(samClient.collectionsForUser(user)) flatMap {
        case Success(collections) =>
          toStrictEntity(Timeout) tflatMap { _ =>
            extractStrictRequest flatMap { request =>
              tprovide((user, collections, request))
            }
          }
        case Failure(e) =>
          throw new RuntimeException(s"Unable to look up collections for user ${user.userId}: ${e.getMessage}", e)
      }
    }
  }

  /**
    * Will verify that none of the GET query parameters are specifying the collection label, and then tack
    * on query parameters for the user's collections on to the query URI
    */
  private def processLabelsForGetQuery(user: User, collections: List[Collection]): Directive1[Uri] = {
    extractUri flatMap { uri =>
      val query = uri.query()

      val labelsWithCollections = processLabels(user, collections, query.getAll(LabelKey))
      val newQueryBuilder = query.newBuilder
      labelsWithCollections foreach { l => newQueryBuilder += (LabelKey -> l) }
      provide(uri.withQuery(newQueryBuilder.result()))
    }
  }

  /**
    * Will verify that none of the POSTed query parameters are specifying the collection label, and then tack
    * on query parameters w/ a label for each of the user's collections. All of the query parameters are rebuilt
    * into a HttpEntity
    */
  private def processLabelsForPostQuery(user: User, collections: List[Collection]): Directive1[HttpEntity.Strict] = {
    import spray.json._
    import DefaultJsonProtocol._

    entity(as[String]) flatMap { paramString =>
      val params = Try(JsonParser(paramString).asInstanceOf[JsArray].elements map { _.asJsObject })
      params match {
        case Success(p) =>
          val (labelParams, otherParams) = p.partition(_.fields.keySet.contains(LabelKey))
          val origLabels = labelParams flatMap { _.fields.values map { _.convertTo[String] } }
          val labels = processLabels(user, collections, origLabels)
          val newLabelParams = labels map { l => JsObject(LabelKey -> JsString(l)) }
          val newParams = JsArray(otherParams ++ newLabelParams)
          provide(HttpEntity(ContentTypes.`application/json`, newParams.toString))
        case Failure(e) => throw InvalidQueryException(e)
      }
    }
  }

  /**
    *
    * Ensures no provided labels include collection names, and then adds the user's collections to the set of labels
    */
  protected[this] def processLabels(user: User, collections: List[Collection], labels: Iterable[String]): Iterable[String] = {
    if (hasCollectionLabel(labels)) {
      log.error("User " + user.userId + " attempted to submit query with label " + CollectionLabelName)
      throw new LabelContainsCollectionException
    } else {
      labels ++ collections.map(c => s"$CollectionLabelName:${c.name}")
    }
  }
}

object QuerySupport {
  def hasCollectionLabel(labels: Iterable[String]): Boolean = labels exists { _.startsWith(s"$CollectionLabelName:") }

  final case class InvalidQueryException(e: Throwable) extends
    Exception(s"Invalid JSON in query POST body: ${e.getMessage}", e)
  val LabelKey = "label"
}
