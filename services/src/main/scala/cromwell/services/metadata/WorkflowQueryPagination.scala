package cromwell.services.metadata

import cromwell.services.metadata.MetadataService.QueryMetadata
import spray.http.HttpHeaders.Link
import spray.http.{HttpHeader, Uri}


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
