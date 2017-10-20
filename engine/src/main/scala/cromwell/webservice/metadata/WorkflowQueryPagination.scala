package cromwell.webservice.metadata

import akka.http.scaladsl.model.Uri.Query
import akka.http.scaladsl.model.headers.{Link, LinkParams}
import akka.http.scaladsl.model.{HttpHeader, Uri}
import cromwell.services.metadata.MetadataService.QueryMetadata


/**
  * Attempts to add query parameters for pagination.
  *
  * NOTE: This is effectively broken, as the returned links are not suitable for use by cromwell clients.
  *
  * This discards the search parameters for GETs, for example it drops parameters such as "start" and "end". Also
  * generates links incompatible with POSTs, as the endpoints read parameters from the HTTP body during POST, __not__
  * from the URI.
  *
  * This may need to receive an entire `HttpRequest` and not just the `Uri` to ensure that it doesn't generate links for POST.
  *
  * The existing `CromwellApiServiceSpec` should be updated to verify the expected behavior for both GET and POST.
  *
  * Left behind for legacy reasons, but don't believe anyone has ever used these non-functional links.
  *
  * Note: As of 6/7/17 the above is confirmed by JG, but leaving it mostly as-is for now
  */
object WorkflowQueryPagination {

  private def generatePaginationParams(page: Int, pageSize: Int): Query = {
    Query(s"page=$page&pagesize=$pageSize")
  }

  //Generates link headers for pagination navigation https://tools.ietf.org/html/rfc5988#page-6
  def generateLinkHeaders(uri: Uri, metadata: Option[QueryMetadata]): List[HttpHeader] = {
    //strip off the query params
    val baseUrl = uri.scheme + ":" + uri.authority + uri.path
    metadata match {
      case Some(meta) =>
        (meta.page, meta.pageSize) match {
          case (Some(p), Some(ps)) =>

            val firstLink = Link(Uri(baseUrl).withQuery(generatePaginationParams(1, ps)), LinkParams.first)

            val prevPage = math.max(p - 1, 1)
            val prevLink = Link(Uri(baseUrl).withQuery(generatePaginationParams(prevPage, ps)), LinkParams.prev)

            val lastPage = math.ceil(meta.totalRecords.getOrElse(1).toDouble / ps.toDouble).toInt
            val lastLink = Link(Uri(baseUrl).withQuery(generatePaginationParams(lastPage, ps)), LinkParams.last)

            val nextPage = math.min(p + 1, lastPage)
            val nextLink = Link(Uri(baseUrl).withQuery(generatePaginationParams(nextPage, ps)), LinkParams.next)

            List(firstLink, prevLink, nextLink, lastLink)
          case _ => List.empty
        }
      case None => List.empty
    }
  }
}
