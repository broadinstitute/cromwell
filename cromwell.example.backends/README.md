# Cromwell Example Backends

This is a folder of example backend providers for Cromwell. You can read about
the providers here, and then copy paste one or more of the providers you want
to use to your Cromwell configuration file, represented here as the
[cromwell.example.conf](../cromwell.example.conf) file in the base of the 
repository.

## What are the backend providers?

### Cloud Providers

 - [AWS](AWS.conf): Amazon Web Services ([documentation](https://cromwell.readthedocs.io/en/stable/tutorials/AwsBatch101/))
 - [BCS](BCS.conf) Alibaba Cloud Batch Compute (BCS) backend ([documentation](https://cromwell.readthedocs.io/en/stable/backends/BCS/))
 - [TES](TES.conf) is a backend that submits jobs to a server with protocol defined by GA4GH ([documentation](https://cromwell.readthedocs.io/en/stable/backends/TES/))
 - [PAPIv2](PAPIv2.conf): Google Pipelines API backend (version 2!) ([documentation](https://cromwell.readthedocs.io/en/stable/backends/Google/))

### Containers

 - [Docker](Docker.conf): an example backend that only runs workflows with docker in *every* command
 - [Singularity](singularity.conf): run Singularity containers locally ([documentation](https://cromwell.readthedocs.io/en/develop/tutorials/Containers/#local-environments))
 - [Singularity+Slurm](singularity.slurm.conf): An example using Singularity with SLURM ([documentation](https://cromwell.readthedocs.io/en/develop/tutorials/Containers/#job-schedulers))
 - [TESK](TESK.conf) is the same, but intended for Kubernetes. See the [TES docs](https://cromwell.readthedocs.io/en/stable/backends/TES/) at the bottom.
 - [udocker](udocker.conf): to interact with udocker locally [documentation](https://cromwell.readthedocs.io/en/develop/tutorials/Containers/#udocker)
 - [udocker+Slurm](udocker.slurm.conf): to interact with udocker on SLURM ([documentation](https://cromwell.readthedocs.io/en/develop/tutorials/Containers/#udocker))

### Workflow Managers

 - [HtCondor](HtCondor.conf): a workload manager at UW-Madison ([documentation](https://cromwell.readthedocs.io/en/stable/backends/HTcondor/))
 - [LSF](LSF.conf): the Platform Load Sharing Facility backend ([documentation](https://cromwell.readthedocs.io/en/stable/backends/LSF/))
 - [SGE](SGE.conf): a backend for Sungrid Engine ([documentation](https://cromwell.readthedocs.io/en/stable/backends/SGE))
 - [slurm](slurm.conf): SLURM workload manager ([documentation](https://cromwell.readthedocs.io/en/stable/backends/SLURM/))
 - [Spark](Spark.conf): a backend for a Spark cluster ([documentation](https://cromwell.readthedocs.io/en/stable/backends/Spark/))

### Custom

 - [LocalExample](LocalExample.conf): What you should use if you want to define a new backend provider ([documentation](https://cromwell.readthedocs.io/en/stable/backends/Local/))


## How do I add a backend provider?

The section in the file called "backends" has a key, "providers" that looks like
this:

```

backend {

  # Override the default backend.
  #default = "LocalExample"

  # The list of providers. Copy paste the contents of a backend provider in this section
  providers {
        ....
  }

}
```

The examples here also have this section. You would want to copy paste the content
of the file, specifically the section for the provider under backend -> providers,
into the backend -> providers section in the [cromwell.example.conf](../cromwell.example.conf).
Here is what it would look like to add the [slurm](slurm.conf) backend
provider example. 

```
backend {

  # Override the default backend.
  #default = "LocalExample"

  # The list of providers. Copy paste the contents of a backend provider in this section
  providers {
    slurm {
         ...
      }
    }

    # Second backend provider would be copy pasted here!

  }
}
```

This isn't json, so you don't need to add commas between the providers - just
copy paste them one after the other in the backend -> providers section.
Let's say we wanted slurm to be our default! We would do this:

```
backend {

  # Override the default backend.
  default = slurm

  # The list of providers. Copy paste the contents of a backend provider in this section
  providers {
    slurm {
         ...
      }
    }
  }
}
```

Don't forget to customize the sections for your purposes! If anything is
not explained clearly, please [open an issue](https://github.com/broadinstitute/cromwell/issues).

## What if a provider is missing?

If a provider is missing and you don't want to use the [LocalExample](LocalExample.conf)
to write a custom provider, please [let us know](https://github.com/broadinstitute/cromwell/issues)
and we can start discussion about how to define your backend.
