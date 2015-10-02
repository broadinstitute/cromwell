package cromwell.util.google

trait ClientSecrets {
  val clientId: String
  val clientSecret: String
}

case class SimpleClientSecrets(clientId: String, clientSecret: String) extends ClientSecrets
