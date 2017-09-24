package cromwell.cloudSupport.gcp.http

import com.google.cloud.http.HttpTransportOptions
import scala.concurrent.duration._

object GoogleHttpTransportOptions {
  val TransportOptions = HttpTransportOptions.newBuilder()
    .setReadTimeout(3.minutes.toMillis.toInt)
    .build()
}
