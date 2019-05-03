<link rel="stylesheet" href="$base$/css/codemirror.css"/>
<link rel="stylesheet" href="$base$/codemirror/addon/display/fullscreen.css"/>
<link rel="stylesheet" href="$base$/codemirror/addon/fold/foldgutter.css" />
<script src="$base$/codemirror/lib/codemirror.js"></script>
<script src="$base$/codemirror/mode/clike/clike.js"></script>
<script src="$base$/codemirror/mode/markdown/markdown.js"></script>
<script src="$base$/codemirror/addon/display/fullscreen.js"></script>
<script src="$base$/codemirror/addon/edit/matchbrackets.js"></script>
<script src="$base$/codemirror/keymap/emacs.js"></script>
<script src="$base$/codemirror/addon/comment/comment.js"></script>
<script src="$base$/codemirror/addon/dialog/dialog.js"></script>
<script src="$base$/codemirror/addon/search/searchcursor.js"></script>
<script src="$base$/codemirror/addon/search/search.js"></script>
<script src="$base$/codemirror/addon/fold/foldcode.js"></script>
<script src="$base$/codemirror/addon/fold/foldgutter.js"></script>
<script src="$base$/codemirror/addon/fold/brace-fold.js"></script>
<script src="$base$/codemirror/addon/fold/xml-fold.js"></script>
<script src="$base$/codemirror/addon/fold/markdown-fold.js"></script>
<script src="$base$/codemirror/addon/fold/comment-fold.js"></script>
<form action="$pageUrl$" method="post" id="editform"
      enctype="application/x-www-form-urlencoded">
<textarea cols="80" name="editedText" id="editedText">
$if(edit)$
$edit$
$endif$
</textarea>
  <br/><label>Description of changes:</label><br/>
  <input type="text" name="logMsg" id="logMsg" value="$logMsg$"/>
  $checkboxes$
  <input type="button" class="editButton" id="update"
  	 onClick="save();" value="Save"/>
  &nbsp;<input type="submit" name="cancel" id="cancel" value="Discard"/>
  &nbsp;<input type="button" class="editButton" id="previewButton"
	       onClick="updatePreviewPane();" style="display: none;"
	       value="Preview"/>
  &nbsp;<input type="button" class="editButton" id="runButton"
	       onClick="run();" style="display: none;"
	       value="Run"/>
  <div id="status"></div>
  <div id="previewpane"></div>
</form>
<script type="text/javascript">
    function setup_codemirror() {
	CodeMirror.keyMap.emacs.Enter = "newlineAndIndent";
	CodeMirror.extendMode("clike", {electricChars: "{};"});
	var editor = 
    	    CodeMirror.fromTextArea(document.getElementById("editedText"), {
		matchBrackets: true,
		mode: "text/x-basiliskc",
		theme: "basilisk",
		tabSize: 8,
		electricChars: true,
		lineWrapping: true,
		foldGutter: true,
		keyMap: "emacs",
      		extraKeys: {
		    "F11": function(cm) {
			cm.setOption("fullScreen", !cm.getOption("fullScreen"));
		    },
		    "Esc": function(cm) {
			if (cm.getOption("fullScreen")) cm.setOption("fullScreen", false);
		    },
		    "Ctrl-E": function(cm){ cm.foldCode(cm.getCursor()); }
		},
		gutters: ["CodeMirror-foldgutter"]
            });
	editor.on("blur", function() {editor.save()});
	update();
	var pending;
	editor.on("change", function() {
	    clearTimeout(pending);
	    pending = setTimeout(update, 400);
	});
	function update() {
	    editor.setOption("mode", /.*\.[ch]/.test("$pageUrl$") ? 
			     "text/x-basiliskc" : "markdown");
	}
    }
</script>
