function checkRunning(url, sec) {
    if ($('#status').html().search("running") != -1) {
	setTimeout (function() {
	    $('#status').load(url, function() { checkRunning (url, Math.min(2*sec, 16000)); });
	}, sec);
    }
    else if (sec > 1000)
	location.reload (true);
}

$(document).ready(function() {
    hide_plots();
    if (location.pathname.endsWith (".c")) {
	var url = location.pathname.concat("?status");
	$('#status').load(url, function(data) { checkRunning (url, 1000); });
    }
});
