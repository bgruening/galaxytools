
var renderTestResults = function(testData) {
	var summary = testData["summary"];
	var numTests = summary["num_tests"];
	var numProblems = summary["num_errors"] + summary["num_failures"] + summary["num_skips"];
	var $overview = $("#overview-content");
	var $progress = $(".progress");
	if(numTests == 0) {
		$overview.addClass("alert").addClass("alert-danger").text("No tests were executed.");
		$progress.append($('<div class="progress-bar progress-bar-warning" role="progressbar" style="width: 100%" />'));
	} else if(numProblems > 0) {
		$overview.addClass("alert").addClass("alert-danger").text("There were problems with " + numProblems + " test(s) out of " + numTests + ".");
		var problemPercent = (numProblems/(1.0 * numTests)) * 100.0;
		var successPercent = 100.0 - problemPercent;
		$progress.append($('<div class="progress-bar progress-bar-success" role="progressbar" style="width: ' + successPercent  +  '%" />'));
		$progress.append($('<div class="progress-bar progress-bar-danger" role="progressbar" style="width: ' + problemPercent  +  '%" />'));
	} else {
		$overview.addClass("alert").addClass("alert-success").text("All " + numTests + " test(s) successfully executed.");
		$progress.append($('<div class="progress-bar progress-bar-success" role="progressbar" style="width: 100%" />'));
	}

	var $sidebar = $("#nav-sidebar-tests");
	for(var index in testData["tests"]) {
		var test = testData["tests"][index];
		var testResult = new TestResult(index, test);
		var rawId = testResult.rawId;

		var panelType = testResult.passed ? "panel-success panel-success-custom" : "panel-danger panel-danger-custom";
		var $panel = $('<div class="panel">');
		$panel.addClass(panelType);

		var $panelHeading = $('<div class="panel-heading">');
		var $panelTitle = $('<div class="panel-title">');
		var $a = $('<a class="collapsed" data-toggle="collapse">');
		$a.attr("id", rawId);
		$a.attr("data-target", "#collapse"  + index);
		var testName = testResult.toolName + " (Test #" + (testResult.testIndex + 1) + (testResult.passed ? "" : ", Failed") + ")";
		$a.text(testName);
		var $navLink = $('<a>').attr('href', '#' + rawId).text(testName)
		if(!testResult.passed) {
			$navLink.addClass("text-danger text-danger-custom");
		} else {
			$navLink.addClass("text-success text-success-custom");
		}
		$sidebar.append($('<li>').append( $navLink ) );
		$panelTitle.append($a)
		$panelHeading.append($panelTitle);

		var $panelBody = $('<div class="panel-body panel-collapse collapse" >');
		$panelBody.attr("id", "collapse" + index);

		var $status = $('<div>').text("status: " + testResult.status);
		$panelBody.append($status);
		if(testResult.problems.length > 0) {
			var $problemsLabel = $('<div>').text("problems: ");
			var $problemsDiv = $('<div style="margin-left:10px;">');
			var $problemsUl = $('<ul>');
			for(var problemIndex in testResult.problems) {
				$problemsUl.append($('<li>').append($('<pre>').text(testResult.problems[problemIndex])));
			}
			$problemsDiv.append($problemsUl);
			$panelBody.append($problemsLabel).append($problemsDiv);
		}
		var $commandLabel = $('<div>command:</div>');
		var $stdoutLabel = $('<div>job standard output:</div>');
		var $stderrLabel = $('<div>job standard error:</div>');
		var $command;
		if(testResult.command !== null) {
			$command = $('<pre class="pre-scrollable" style="margin-left:10px;">').text(testResult.command);
		} else {
			$command = $('<div class="alert alert-warning" style="margin-left:10px;">').text("No command recorded.");
		}
		var $stdout;
		if(testResult.stdout !== null) {
			$stdout = $('<pre class="pre-scrollable" style="margin-left:10px;">').text(testResult.stdout);
		} else {
			$stdout = $('<div class="alert alert-warning" style="margin-left:10px;">').text("No standard output recorded.");
		}
		var $stderr;
		if(testResult.stderr !== null) {
			$stderr = $('<pre class="pre-scrollable" style="margin-left:10px;">').text(testResult.stderr);
		} else {
			$stderr = $('<div class="alert alert-warning" style="margin-left:10px;">').text("No standard error recorded.");
		}
		$panelBody
			.append($commandLabel)
			.append($command)
			.append($stdoutLabel)
			.append($stdout)
			.append($stderrLabel)
			.append($stderr);
		if(!testResult.passed) {
			var $logLabel = $('<div>log:</div>');
			var $log = $('<pre class="pre-scrollable" style="margin-left: 10px;">').text(testResult.problemLog);
			$panelBody.append($logLabel).append($log);
		}

		$panel.append($panelHeading).append($panelBody);
		$(".main").append($panel);
	}
}

var TestResult = function(index, data) {
	this.rawId = data["id"];

	var idParts = this.rawId.split("TestForTool_");
	var testMethod = idParts[idParts.length-1];
	var splitParts;
	if(testMethod.indexOf(".test_tool_") > -1) {
		splitParts = testMethod.split(".test_tool_");
	} else {
		splitParts = rSplit(testMethod, "-", 1);
	}
	var toolName = splitParts[0];
	var testIndex;
	if(data["data"]["test_index"] !== null) {
		testIndex = data["data"]["test_index"];
	} else {
		testIndex = splitParts[1];
	}
	this.toolName = toolName;
	this.testIndex = parseInt(testIndex === undefined ? index : testIndex);
	this.status = data["data"]["status"];
	var job = data["data"]["job"];
	if(job) {
		this.stdout = data["data"]["job"]["stdout"];
		this.stderr = data["data"]["job"]["stderr"];
		this.command = data["data"]["job"]["command_line"];
	} else {
		this.stdout = null;
		this.stderr = null;
		this.command = null;
	}
	this.problems = [];
	var outputProblems = data["data"]["output_problems"] || [];
	var executionProblem = data["data"]["execution_problem"];
	this.problems.push.apply(this.problems, outputProblems);
	if(executionProblem) {
		this.problems.push(executionProblem);
	}
	this.problemLog = data["data"]["problem_log"];
	this.passed = (this.status == "success");
}

// http://stackoverflow.com/questions/5202085/javascript-equivalent-of-pythons-rsplit
function rSplit(str, sep, maxsplit) {
    var split = str.split(sep);
    return maxsplit ? [ split.slice(0, -maxsplit).join(sep) ].concat(split.slice(-maxsplit)) : split;
}


// http://stackoverflow.com/questions/19491336/get-url-parameter-jquery
function getUrlParameter(sParam)
{
    var sPageURL = window.location.search.substring(1);
    var sURLVariables = sPageURL.split('&');
    for (var i = 0; i < sURLVariables.length; i++)
    {
        var sParameterName = sURLVariables[i].split('=');
        if (sParameterName[0] == sParam)
        {
            return sParameterName[1];
        }
    }
}
