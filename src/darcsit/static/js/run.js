function run() {
    if (!check_path())
	return;

    $('#status').html("<br><div class=message><div id=msg_logo><img src=/img/gears.gif></div><div id=msg_label>The code is waiting to run.</div></div>");
    $('#runButton').attr('disabled','disabled');

    var url = location.pathname.replace(/_edit\//,"_run/");
    $.post(
	url,
	$("#editedText").serialize(),
	function(data) {
	    $('#messages').html(data);
	    updatePreviewPane();
	    $('#runButton').removeAttr('disabled');
	},
	"html");
};

$(document).ready(function(){
    $("#runButton").show();
});
