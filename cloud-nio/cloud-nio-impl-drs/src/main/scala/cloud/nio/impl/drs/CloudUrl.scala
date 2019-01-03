package cloud.nio.impl.drs



sealed trait CloudUrl

case class GcsUrl(url: String) extends CloudUrl {
}
