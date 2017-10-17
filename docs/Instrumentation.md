## StatsD

Cromwell collects metrics while it's running and sends them to an internal service. By default this service ignores those metrics, but it can be configured to forward them to a [StatsD](https://github.com/etsy/statsd) server.
To do so, simply add this snippet to your configuration file:

```hocon
services.Instrumentation.class = "cromwell.services.instrumentation.impl.akka.AkkaInstrumentationServiceActor"
```
Make sure to configure your StatsD service:

```hocon
services.Instrumentation.config.statsd {
    hostname = "localhost" # Replace with your host
    port = 8125 # Replace with your port
    # prefix = "my_prefix" # All metrics will be prefixed by this value if present.
    flush-rate = 1 second # Rate at which metrics are sent to the StatsD server
  }
```

There is also an additional configuration value that can be set: 

```hocon
# Rate at which Cromwell updates its gauge values (number of workflows running, queued, etc...)
system.instrumentation-rate = 5 seconds
```

### Metrics

The current StatsD implementation uses metrics-statsd to report instrumentation values.
metrics-statsd reports all metrics with a gauge type.
This means all the metric will be under the gauge section. We might add /remove metrics in the future depending on need and usage.
However here the current high level categories:

`backend`, `rest-api`, `job`, `workflow`, `io`