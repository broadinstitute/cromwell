## Configuration Files

### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Five Minute Introduction](FiveMinuteIntro.md)


### Goals

At the end of this tutorial you'll have set up a configuration file for Cromwell and used it to modify Cromwell's behavior.

### Let's get started

#### Customizing Cromwell with Configuration Files

When Cromwell runs, it contains a large number of default options useful for getting started. For example, by default Cromwell doesn't require an external database while running all workflow jobs on your local machine.

Soon you may want to start storing the results of your Cromwell runs in an external MySQL database. Or, you may want to run jobs on your organizations compute farm, or even run jobs in the cloud via the Pipelines API. All of these changes to the defaults will done by setting configuration values.

When you have many configuration settings you would like to set, you specify them in custom configuration file. See the [configuration](../Configuring) page for more specific information on the configuration file, and for links to the example configuration file.

#### Configuration file syntax

Cromwell configuration files are written in a syntax called HOCON. See the [HOCON documentation](https://github.com/typesafehub/config/blob/master/HOCON.md#hocon-human-optimized-config-object-notation) for more information on all the ways one can create a valid configuration file.

#### Creating your first configuration file

To get started customizing Cromwell via a configuration file, create a new empty text file, say `your.conf`. Then add this include at the top:

```hocon
include required(classpath("application"))
```

The default Cromwell configuration values are set via Cromwell's `application.conf`. To ensure that you always have the defaults from the `application.conf`, you must include it at the top of your new configuration file.

#### Running Cromwell with your configuration file

One you have created a new configuration file, you can pass the path to Cromwell by setting the system property `config.file`:

```bash
java -Dconfig.file=/path/to/your.conf -jar cromwell-[VERSION].jar server
```

Cromwell should start up as normal. As you haven't actually overridden any values yet, Cromwell should be running with the same settings.

#### Setting a configuration value

To override a configuration value, you can specify new values in your configuration file. For example, say you want to change the default port that cromwell listens from `8000` to `8080`. In your config file you can set:

```hocon
# below the include line from before
webservice {
  port = 8080
}
```

When you then run Cromwell updated config file, cromwell will now be listening on 8080 or 8000.

#### Finding more configuration properties

In addition to the common configuration properties listed on the [configuration](../Configuring) page, there are also a large number of example configuration stanzas commented in [cromwell.examples.conf][cromwell-examples-conf], and
backend provider examples in [cromwell.example.backends][cromwell-examples-folder].


[cromwell.examples.conf](https://github.com/broadinstitute/cromwell/blob/develop/cromwell.examples.conf).

### Next Steps

After completing this tutorial you might find the following page interesting:

* [Configuring the Local Backend](LocalBackendIntro)
* [Server Mode](ServerMode.md)

[cromwell-examples-conf]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.examples.conf
[cromwell-examples-folder]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends
