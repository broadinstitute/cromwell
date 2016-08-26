package cromwell.webservice.metadata

import cromwell.services.metadata.MetadataService.QueryMetadata
import spray.http.HttpHeaders.Link
import spray.http.{HttpHeader, Uri}

/**
  * Attempts to add query parameters for pagination.
  *
  * NOTE: This trait is effectively broken, as the returned links are not suitable for use by cromwell clients.
  *
  * The trait discards the search parameters for GETs, for example it drops parameters such as "start" and "end". Also
  * generates links incompatible with POSTs, as the endpoints read parameters from the HTTP body during POST, __not__
  * from the URI.
  *
  * This trait may need to receive an entire `spray.http.HttpRequest` and not just the `spray.http.Uri` to ensure that
  * it doesn't generate links for POST.
  *
  * The existing `CromwellApiServiceSpec` should be updated to verify the expected behavior for both GET and POST.
  *
  * Left behind for legacy reasons, but don't believe anyone has ever used these non-functional links.
  */
trait WorkflowQueryPagination {

  protected def generatePaginationParams(page: Int, pageSize: Int): String = {
    s"page=$page&pagesize=$pageSize"
  }

  //Generates link headers for pagination navigation https://tools.ietf.org/html/rfc5988#page-6
  protected def generateLinkHeaders(uri: Uri, metadata: Option[QueryMetadata]): Seq[HttpHeader] = {
    //strip off the query params
    val baseUrl = uri.scheme + ":" + uri.authority + uri.path
    metadata match {
      case Some(meta) =>
        (meta.page, meta.pageSize) match {
          case (Some(p), Some(ps)) =>

            val firstLink = Link(Uri(baseUrl).withQuery(generatePaginationParams(1, ps)), Link.first)

            val prevPage = math.max(p - 1, 1)
            val prevLink = Link(Uri(baseUrl).withQuery(generatePaginationParams(prevPage, ps)), Link.prev)

            val lastPage = math.ceil(meta.totalRecords.getOrElse(1).toDouble / ps.toDouble).toInt
            val lastLink = Link(Uri(baseUrl).withQuery(generatePaginationParams(lastPage, ps)), Link.last)

            val nextPage = math.min(p + 1, lastPage)
            val nextLink = Link(Uri(baseUrl).withQuery(generatePaginationParams(nextPage, ps)), Link.next)

            Seq(firstLink, prevLink, nextLink, lastLink)

          case _ => Seq()
        }
      case None => Seq()
    }
  }
}
