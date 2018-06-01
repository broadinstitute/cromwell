Workflow options files rendered from Vault Consul Templates (.ctmpl files) will be copied to this directory by
Cromwell's CI framework prior to building the Centaur jar. For those Centaur builds that don't actually do ctmpl
rendering, dummy versions of the options files should be committed here so Centaur won't fail to run if it can't find
required workflow options files.
