function hide_plots()
{
    $("[id^=plot]").each(function( index, element ){
	id = $(this).attr('id');
	$('#after' + id).append( $(this) );
	$('#after' + id).hide();
	$('#button' + id).click( function(e) {
	    e.preventDefault();
	    $($(this).attr('id').replace('button', '#after')).toggle();
	    return false;
	});
    });
}
