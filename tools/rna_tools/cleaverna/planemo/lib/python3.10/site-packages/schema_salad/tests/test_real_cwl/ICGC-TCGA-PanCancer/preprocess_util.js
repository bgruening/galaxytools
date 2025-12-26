
function filterForIndels(inArr)
{
	var arr = [];
	for (var i = 0; i < inArr.length ; i++)
	{
		if (inArr[i].basename.indexOf("indel") >= 0)
		{
			arr.push(inArr[i]);
		}
	}
	return arr;
}

/**
 * workflowName - the name of the workflow to filter for.
 * vcfType - the type of VCF (snv, indel, etc...)
 * inArr - an array of files (File[]) to search through.
 */
function filterFor(workflowName, vcfType, inArr)
{
	var arr = [];
	for (var i = 0; i < inArr.length; i++)
	{
		if (typeof(inArr[i]) == "string") //("class" in inArr[i] && inArr[i].class == "File")
		{
			if (inArr[i].indexOf(workflowName) >= 0 && inArr[i].indexOf(vcfType) >= 0)
			{
				arr.push(inArr[i])
			}
		}
		else
		{
			if ("class" in inArr[i] && inArr[i].class == "File")
			{
				if (inArr[i].basename.indexOf(workflowName) >= 0 && inArr[i].basename.indexOf(vcfType) >= 0)
				{
					arr.push(inArr[i])
				}
			}
		}
	}
	return arr;
}
