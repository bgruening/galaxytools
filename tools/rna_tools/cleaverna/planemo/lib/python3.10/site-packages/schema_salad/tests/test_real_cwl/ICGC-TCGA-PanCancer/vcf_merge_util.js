
function formatArray(arr)
{
	var files = "";
	for (var i = 0; i < arr.length; i++)
	{
		files = files.concat(arr[i].path,",");
	}
	return files;
}
