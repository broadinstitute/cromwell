# Overview

Cromwell's instrumentation support can be useful to collect utilization data in long-running, high-volume
production environments. The default implementation of this ignores these metrics, but Cromwell includes alternate implementations that can forward metrics to a
specific server.

### StatsD

While this instrumentation support can be used in smaller environments it will still require setting up a
[StatsD](https://github.com/etsy/statsd) server outside of Cromwell and it's possible not enough data would be produced to be useful. 
Cromwell collects metrics while running and sends them to an internal service. 

Make sure to configure your StatsD service:

```hocon
services.Instrumentation {
	class = "cromwell.services.instrumentation.impl.statsd.StatsDInstrumentationServiceActor"

	config {
	    hostname = "localhost" # Replace with your host
	    port = 8125 # Replace with your port
	    # prefix = "my_prefix" # All metrics will be prefixed by this value if present.
	    flush-rate = 1 second # Rate at which metrics are sent to the StatsD server
  	}
}
```

There is also an additional configuration value that can be set: 

```hocon
# Rate at which Cromwell updates its gauge values (number of workflows running, queued, etc...)
system.instrumentation-rate = 5 seconds
```

If you have multiple Cromwell instances, and would like to separate the instrumentation path for each instance, set the `system.cromwell_id` with the unique identifier for each Cromwell instance. For example,
```hocon
system.cromwell_id = "cromwell-instance-1"
```
will prepend all the metrics with path `cromwell.cromwell-instance-1...` for that instance.


##### Metrics

The current StatsD implementation uses metrics-statsd to report instrumentation values.
metrics-statsd reports all metrics with a gauge type.
This means all metrics will be under the gauge section. We might add or remove metrics in the future depending on need and usage.
These are the current high level categories:

* `backend`
* `rest-api`
* `job`
* `workflow`
* `io`


### Stackdriver

Cromwell now supports sending metrics to [Google's Stackdriver API](https://cloud.google.com/monitoring/api/v3/). To use the Stackdriver instrumentation
specify this in your config:
```hocon
services.Instrumentation {
    class = "cromwell.services.instrumentation.impl.stackdriver.StackdriverInstrumentationServiceActor"

    config {
        # auth scheme can be `application_default` or `service_account`
        auth = "service-account"
        google-project = "my-project"
        # rate at which aggregated metrics will be sent to Stackdriver. It needs to be equal or greater than 1 minute.
        # Google's Stackdriver API needs each metric to be sent not more than once per minute.
        flush-rate = 1 minute
        # below 3 keys are attached as labels to each metric. `cromwell-perf-test-case` is specifically meant for perf env.
        cromwell-instance-role = "role"
        cromwell-perf-test-case = "perf-test-1"
    }
 }
```
The 2 label keys are optional. If specified, each metric will have label(s) added in the form of a (key, value) pair.
So for example, if `cromwell-instance-role = "backend"` is mentioned in config, each metric data point sent to Stackdriver
will have a label (cromwell_instance_role, backend) added to it.

There is another optional label that can be added to each metric. `cromwell_id` represents the identifier for different Cromwell instances.
```hocon
# Unique Cromwell instance identifier
system.cromwell_id = "cromwell-instance-1"
```


##### Metric type and Label keys naming convention
More details on the this can be found [here](https://cloud.google.com/monitoring/api/v3/metrics-details#metric-kinds).

You must adhere to the following spelling rules for metric type names:
- You can use upper and lower-case letters, digits, and underscores (_) in the names.
- You can use periods (.) in the domain part of the names.
- You can use forward slashes (/) to separate path elements.
- You can start each path element with a letter or digit.
- The maximum length of a metric type name is 200 characters.

You must adhere to the following spelling rules for metric label names:
- You can use upper and lower-case letters, digits, underscores (_) in the names.
- You can start names with a letter or digit.
- The maximum length of a metric label name is 100 characters.



