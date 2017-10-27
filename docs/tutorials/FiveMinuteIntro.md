# Five minute Introduction to Cromwell

### Prerequisites:

* A Unix-based operating system (yes, that includes Mac!)
* A Java 8 runtime environment 
	* You can see what you have by running `java --version` on a terminal. You're looking for a version that's at least `1.8` or higher.
	* If not, you can download Java [here](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)
* A sense of adventure!

### Goals

At the end of this five minute introduction you will have:

- Downloaded Cromwell!
- Written your first workflow
- Run it through Cromwell

### Step 1: Downloading Cromwell

We host our Cromwell releases on GitHub! You can find the latest version on our [Releases](https://github.com/broadinstitute/cromwell/releases/latest) page.

* Head down to the bottom of the page and find the link to download `cromwell-29.jar`. Note that when you read these instructions, the latest Cromwell version might be higher than 29. That's fine, these instructions will work just as well! Download the jar file.

* You'll probably want to put your downloaded Cromwell somewhere you can find it later. If you want, you can make a cromwell directory in your home directory and put it there.

* For example, in a terminal (if you're not using a Mac, the final command might be different for you!):
```sh
$ cd ~
$ mkdir cromwell
$ cd cromwell
$ cp ~/Downloads/cromwell-29.jar .
```

### Step 2: Writing your first workflow description

This bit's easy, you're just going to copy and paste something from the internet!

Create a file in your new `cromwell` directory called `myWorkflow.wdl`. Copy the following content into it.

```wdl
workflow myWorkflow {
	call myTask
}

task myTask {
	command {
		echo "hello world"
	}
	output {
		String out = read_string(stdout())
	}
}
```

Don't worry too much about the workflow contents for now. If you'd like to learn more about how to author WDL, you can hopefully find all the resources you could ever want [here](https://github.com/openwdl/wdl)!


### Next Steps

After completing this tutorial you might find the following page interesting:

* [ServerMode](ServerMode)
* [Configuration Files](ConfigurationFiles)
