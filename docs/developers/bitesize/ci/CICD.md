Cromwell uses [Continuous Integration](https://en.wikipedia.org/wiki/Continuous_integration) (CI) testing, along with [Continuous Delivery](https://en.wikipedia.org/wiki/Continuous_delivery) (CD) to the Cromwell-as-a-Service (CaaS) `DEV` environment. [Continuous Deployment](https://en.wikipedia.org/wiki/Continuous_deployment) is not implemented.

## CI Testing in Github Actions

All CI testing of Cromwell is done through [Github Actions](https://github.com/broadinstitute/cromwell/actions). Github Actions will test every pull request submitted by a trusted contributor. Pull requests must pass all unit and integration tests in order to be merged into the develop environment. 

## Vulnerability Scanning

In collaboration with the DSP Information Security team, scans include but are not limited to:

- Committed Git secrets
- Vulnerable Java dependencies
- Penetration testing
