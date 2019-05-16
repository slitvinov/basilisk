<div id="content">
  $if(revision)$
    <h2 class="revision">Revision $revision$</h2>
  $endif$
  $if(notitle)$
  $else$
    <h1 class="pageTitle"><a href="$pageUrl$">$pagetitle$</a></h1>
  $endif$
  $messages()$
  $if(editable)$
    $edit()$
  $else$
    <div id="status"></div>
  $endif$
  $if(login)$
    $login()$
  $endif$
  $if(register)$
    $register()$
  $endif$
  $content$
</div>