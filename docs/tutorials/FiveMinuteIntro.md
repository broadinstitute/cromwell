# Five minute Introduction to Cromwell

### Prerequisites:

* A Unix-based operating system (yes, that includes Mac!)
* A Java 8 runtime environment 
	* You can see what you have by running `$ java -version` on a terminal. You're looking for a version that's at least `1.8` or higher.
	* If not, you can download Java [here](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html).
* A sense of adventure!

### Goals

At the end of this five minute introduction you will have:

- Downloaded Cromwell!
- Written your first workflow
- Run it through Cromwell

### Step 1: Downloading Cromwell

We host our Cromwell releases on GitHub! You can find the latest version on our [Releases](https://github.com/broadinstitute/cromwell/releases/latest) page.

* Look for the latest version at the top of the page, and find the link to download the jar. It'll have a name like `cromwell-XY.jar` where `XY` is the version. Download the jar file.
* WARNING! If you're on a Mac, the security settings might try to stop you from running Cromwell! Don't worry, if this happens just go to `System Preferences > Security & Privacy > General` and find the `cromwell` jar listed on the page. Click `Open anyway`. The `cromwell-XY.jar` will now automatically download to your `Downloads` directory.
* Put your downloaded Cromwell somewhere you can find it later, like in a Cromwell directory in your home directory.

For example, in a terminal:
```sh
cd ~
mkdir cromwell
cp ~/Downloads/cromwell-XY.jar cromwell/
cd cromwell/
```
_(if you're not using a Mac, the final command might be different for you)_


### Step 2: Writing your first workflow description

This bit is easy, you're just going to copy and paste something from the internet.

Open your favorite editor. Paste in the following content and save it as `myWorkflow.wdl` in your new `cromwell` directory:

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

Don't worry, **you don't need to understand too much about the workflow contents to continue for now**. In brief, it tells Cromwell to run a task to run `echo "hello world"`, and then return the output as a String. If you'd like to learn more about how to author WDL, you can find all the WDL resources you could ever want [here](https://github.com/openwdl/wdl).

### Step 3: Running the workflow

Ok, we have Cromwell, we have a workflow, let's put it all together! 

Make sure you're in the cromwell directory with the `.jar` file and the `.wdl` file. Now type in:
```sh
java -jar cromwell-XY.jar run myWorkflow.wdl
```

Cromwell will print out a fair old chunk of logging information, which can be configured (once you've completed this tutorial and [Configuration Files](ConfigurationFiles), you might want to investigate the [Logging](../Logging) page)

Ultimately, the workflow should succeed and you'll end up with the following output printed out when Cromwell finishes:
```json
{
	"myWorkflow.myTask.out": "hello world"
}
```

Ok, you can stop your timer! You just installed and ran your first workflow in Cromwell, congratulations!

### Next Steps

Pat yourself on the back for completing this tutorial, bravo! Then continue on to one of the follow pages:

* [Server Mode](ServerMode)
* [Configuration Files](ConfigurationFiles)
