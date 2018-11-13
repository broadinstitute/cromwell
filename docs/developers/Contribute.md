**Contribution Guidelines**

Cromwell is an open-source project and we heartily welcome community contributions to both our code and our documentation. Here are some guidelines for adding documentation and recommendations on where we could use help the most. 

First off, here are useful links:

* Cromwell documentation: [cromwell.readthedocs.io](http://cromwell.readthedocs.io)
* Source on Github: [github.com/broadinstitute/cromwell](https://github.com/broadinstitute/cromwell/tree/develop/docs)
* Builds on ReadTheDocs: [readthedocs.org/projects/cromwell](https://readthedocs.org/projects/cromwell/builds/)
* How to build and view the documentation locally: [mkdocs.readthedocs.io](https://mkdocs.readthedocs.io/en/stable/#installation)
* Public Google Docs - Table of Contents:  [drive.google.com/open?id=1myTVzWx5HG720nPBUAF9vE8XTacrzCW70RpgsjMnchM](https://drive.google.com/open?id=1myTVzWx5HG720nPBUAF9vE8XTacrzCW70RpgsjMnchM)

### "needs docs"

There are plenty of areas of the Cromwell documentation that need to be updated, improved, or added. Within the Cromwell repo on Github there are many issues labeled as ["needs docs"](https://github.com/broadinstitute/cromwell/issues?q=is%3Aopen+is%3Aissue+label%3A%22needs+docs%22), so feel free to start there. 

If you would like to request additional documentation, you can create a [Github issue in the Cromwell repo](https://github.com/broadinstitute/cromwell/issues/new).

### Writing Tips

1. Keep it clear, accurate, and concise.
2. Put the most important information first.
3. Use the second person, use “you” instead of “the user”.
4. No passive verbs (everything is done by something).
5. Link to the original source, don't repeat documentation. 

### Formatting

The documentation is written in Markdown. Click here for a [Github Guide on Markdown](https://guides.github.com/features/mastering-markdown/), and click here for more [tips from MkDocs](http://www.mkdocs.org/user-guide/writing-your-docs/).

### Styling

**Links:**

* Absolute: `[link text](www.destinationURL.com)`  
_Example:_ `[Broad Institute](www.broadinstitute.org)` _produces this link_ [Broad Institute](www.broadinstitute.org).
* Relative: `[link text](Destination_Page)`, where `Destination_Page` is the file name without the `.md` extension  
_Example:_ `[How to use the Cromwell CLI](CommandLine)` _produces this link_ [How to use the Cromwell CLI](CommandLine).
* Anchor link: `[anchor text](../Path/To/Page#Anchor)`  
_Example:_ `[HPC filesystems](backends/HPC#filesystems)` _produces this link_ [HPC filesystems](backends/HPC#filesystems).

**Code:**

* To style a word of code, use a backtick (\`) before and after the word.  
_Example:_ \`file.json\` _produces_ `file.json`. 
* To style a block of code, use three backticks (\`\`\`) before and after the block of code.  
_Example:_  
	\`\`\`  
		workflow myWorkflow {  
		call myTask  
		}  
	\`\`\`  
_produces this block_  
```  
	workflow myWorkflow {  
	call myTask  
	}  
```  
* To use syntax highlighting, include the language after the first three backticks (\`\`\`).  
_Example:_   
	\`\`\`json  
		\{  
		    "MyWorkflow.MyTask.VariableTwo": "Variable2"  
		\}  
	\`\`\`	
_produces this block_  
```json  
	{
	    "MyWorkflow.MyTask.VariableTwo": "Variable2"
	}
```	

**Images**

* Relative: `![](ImgName.png)`
	* _Example:_ `![](jamie_the_cromwell_pig.png)` _produces this image_  
	![](jamie_the_cromwell_pig.png)  
* Absolute: `![](URLofImg.png)`
	* _Example:_ `![](https://www.broadinstitute.org/sites/all/themes/custom/at_broad/logo.png)` _produces this image_  
	![](https://www.broadinstitute.org/sites/all/themes/custom/at_broad/logo.png)  

### REST API & Menu

**REST API:**

1. Edit the [`cromwell.yaml`](https://github.com/broadinstitute/cromwell/blob/develop/engine/src/main/resources/swagger/cromwell.yaml) to make any changes to the REST API content.  
2. Regenerate the REST API markdown file by running `sbt generateRestApiDocs` from the main Cromwell directory.
3. Commit both the changes to the `cromwell.yaml` and the (re)generated [RESTAPI.md](https://github.com/broadinstitute/cromwell/blob/develop/docs/api/RESTAPI.md).
4. Once your branch is merged to the [`develop` branch](https://github.com/broadinstitute/cromwell/tree/develop), you will see your changes on the [REST API page](api/RESTAPI/).

**Left-side menu:**

To add or remove items from the menu, edit [mkdocs.yml](https://github.com/broadinstitute/cromwell/blob/develop/mkdocs.yml) in Cromwell.

### FAQs

**_Why isn't my documentation showing up?_** 

* **Is your PR merged?**  
If not, [kindly ask the team to merge it](https://github.com/broadinstitute/cromwell/pulls). Once your PR is merged to develop, it will trigger an automatic build.  

* **Has the build finished?**   
[Check build status here](https://readthedocs.org/projects/cromwell/builds/).

* **Did you add the file(s) to the YAML file?**  
If not, [add it here](https://github.com/broadinstitute/cromwell/blob/develop/mkdocs.yml).

**_How do I add to the REST API documentation?_**

* **Don't forget to regenerate**  
After you edit the [`cromwell.yaml`](https://github.com/broadinstitute/cromwell/blob/develop/engine/src/main/resources/swagger/cromwell.yaml), run `sbt generateRestApiDocs` and commit all changes.  
_Hint:_ Once you have regenerated the docs correctly, the hidden timestamp at the top of the [`RESTAPI.md` file](https://raw.githubusercontent.com/broadinstitute/cromwell/develop/docs/api/RESTAPI.md) will show the current time.
