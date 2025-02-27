**Contribution Guidelines**

Cromwell is an open-source project and we heartily welcome community contributions to both our code and our documentation. Here are some guidelines for adding documentation and recommendations on where we could use help the most. 

First off, here are useful links:

* Cromwell documentation: [cromwell.readthedocs.io](http://cromwell.readthedocs.io)
* Source on Github: [github.com/broadinstitute/cromwell](https://github.com/broadinstitute/cromwell/tree/develop/docs)
* Builds on ReadTheDocs: [readthedocs.org/projects/cromwell](https://readthedocs.org/projects/cromwell/builds/)
* How to build and view the documentation locally: [mkdocs.readthedocs.io](https://mkdocs.readthedocs.io/en/stable/#installation)

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
_Example:_ `[Broad Institute](www.broadinstitute.org)` _produces this link_ [Broad Institute](https://www.broadinstitute.org).
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
	* _Example:_ `![](../jamie_the_cromwell_pig.png)` _produces this image_  
	![](../jamie_the_cromwell_pig.png)  
* Absolute: `![](URLofImg.png)`
	* _Example:_ `![](https://www.broadinstitute.org/sites/all/themes/custom/at_broad/logo.png)` _produces this image_  
	![](https://www.broadinstitute.org/sites/all/themes/custom/at_broad/logo.png)  

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
