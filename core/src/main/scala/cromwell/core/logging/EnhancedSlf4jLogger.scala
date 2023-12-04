package cromwell.core.logging

import akka.event.slf4j.Slf4jLogger

class EnhancedSlf4jLogger extends Slf4jLogger {

  /**
    * Format the timestamp as a simple long. Allows the akkaTimestamp to be retrieved later from the MDC by custom
    * converters.
    *
    * NOTE: Should not be necessary once this issue is resolved: 
    *   - https://github.com/akka/akka/issues/18079#issuecomment-125175884
    *
    * @see [[EnhancedDateConverter.convert()]]
    */
  override protected def formatTimestamp(timestamp: Long): String = String.valueOf(timestamp)
}
