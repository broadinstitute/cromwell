[![codecov](https://codecov.io/gh/broadinstitute/cromwell/branch/develop/graph/badge.svg)](https://codecov.io/gh/broadinstitute/cromwell)

## Welcome to Cromwell

Cromwell is an open-source Workflow Management System for bioinformatics. Licensing is [BSD 3-Clause](LICENSE.txt).

The [Cromwell documentation has a dedicated site](https://cromwell.readthedocs.io/en/stable).

First time to Cromwell? Get started with [Tutorials](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/).

### Community

Thinking about contributing to Cromwell? Get started by reading our [Contributor Guide](CONTRIBUTING.md).

Cromwell has a growing ecosystem of community-backed projects to make your experience even better! Check out our [Ecosystem](https://cromwell.readthedocs.io/en/stable/Ecosystem/) page to learn more.

Talk to us:
- [Join the Cromwell Slack workspace](https://join.slack.com/t/cromwellhq/shared_invite/zt-dxmmrtye-JHxwKE53rfKE_ZWdOHIB4g) to discuss the Cromwell workflow engine.
- [Join the OpenWDL Slack workspace](https://join.slack.com/t/openwdl/shared_invite/zt-ctmj4mhf-cFBNxIiZYs6SY9HgM9UAVw) to discuss the evolution of the WDL language itself.
    - More information about WDL is available in [that project's repository](https://github.com/openwdl/wdl).  

### Capabilities and roadmap

Many users today run their WDL workflows in [Terra](https://support.terra.bio/hc/en-us/articles/360036379771-Get-started-running-workflows), a managed bioinformatics platform with built-in Cromwell support.

Users with specialized needs who wish to install and maintain their own Cromwell instances can [download](https://github.com/broadinstitute/cromwell/releases) a JAR or Docker image. The development team accepts reproducible bug reports from self-managed instances, but cannot provide direct support.

Cromwell uses [modular backends](https://cromwell.readthedocs.io/en/stable/backends/Backends/) to support different cloud vendors. The team is currently developing for AWS Batch and GCP Batch. Maintenance of other backends is community-based.

Cromwell [supports](https://cromwell.readthedocs.io/en/stable/LanguageSupport/) the WDL workflow language. Cromwell version 80 and above no longer support CWL.

### Security reports

If you believe you have found a security issue please contact `infosec@broadinstitute.org`.

### Issue tracking

Need to file an issue? Head over to [Github Issues](https://github.com/broadinstitute/cromwell/issues).

![Jamie, the Cromwell pig](docs/jamie_the_cromwell_pig.png)
