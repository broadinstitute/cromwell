package cromwell.instrumentation

import com.typesafe.config.ConfigFactory
import lenthall.config.ScalaConfig._
import kamon.Kamon
import kamon.metric.instrument.Histogram.Record
import kamon.metric.instrument._
import kamon.spray.KamonTraceDirectives
import spray.routing._

/**
 * Bit of a hack. It's not always easy to simply wrap kamon calls with an if/else (and even if it were, it'd get
 * messy quickly). Instead we're setting up a single global monitoring object which is either a simple wrapper around
 * kamon or a noop construct.
 *
 * Note that not all of kamon is wrapped here, just the bits of instrumentation that we're using. If one adds new
 * instrumentation it might require adding stuff here
 */
object Instrumentation {
  private val Config = ConfigFactory.load.getConfigOption("instrumentation")
  private val UseKamon = Config flatMap { _.getBooleanOption("use-kamon") } getOrElse false

  val Monitor = if (UseKamon) KamonInstrumentation else NoOpInstrumentation
}

sealed trait Instrumentation {
  def minMaxCounter(name: String): MinMaxCounter
  def histogram(name: String): Histogram
  def traceName(name: String): Directive0
  def start(): Unit
}

case object KamonInstrumentation extends Instrumentation {
  override def traceName(name: String): Directive0 = KamonTraceDirectives.traceName(name)
  override def histogram(name: String): Histogram = Kamon.metrics.histogram(name)
  override def minMaxCounter(name: String): MinMaxCounter = Kamon.metrics.minMaxCounter(name)
  override def start() = Kamon.start()
}

case object NoOpInstrumentation extends Instrumentation {
  override def traceName(name: String): Directive0 = KamonTraceDirectives.pass
  override def histogram(name: String): Histogram = NoOpHistogram
  override def minMaxCounter(name: String): MinMaxCounter = NoOpMinMaxCounter
  override def start() = ()

  case object NoOpMinMaxCounter extends MinMaxCounter {
    override def increment() = ()
    override def increment(times: Long) = ()
    override def decrement() = ()
    override def decrement(times: Long) = ()
    override val refreshValues = ()
    override val cleanup = ()
    override def collect(context: CollectionContext) = new NoOpInstrumentHistogramSnapshot
  }

  case object NoOpHistogram extends Histogram {
    override def record(value: Long) = ()
    override def record(value: Long, count: Long) = ()
    override val cleanup = ()
    override def collect(context: CollectionContext) = new NoOpInstrumentHistogramSnapshot
  }

  class NoOpInstrumentHistogramSnapshot extends Histogram.Snapshot {
    override val isEmpty = true
    override val numberOfMeasurements = 0L
    override val min = 0L
    override val max = 0L
    override val sum = 0L
    override def percentile(percentile: Double) = 0L
    override val recordsIterator: Iterator[Record] = Iterator.empty
    override def merge(that: InstrumentSnapshot, context: CollectionContext): InstrumentSnapshot = that
    override def merge(that: Histogram.Snapshot, context: CollectionContext): Histogram.Snapshot = that
  }
}


