package cromwell

package object webservice {

  case class QueryParameter(key: String, value: String)
  type QueryParameters = Seq[QueryParameter]

}
