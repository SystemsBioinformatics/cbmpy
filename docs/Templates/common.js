function emcode(name, htext, optdom)
  //Based on SUN IT's script but modified to work more smoothly. Email address now displayed in status bar when
  //onMouseOver without needing to click on url - brett 20070623
  {
    var site='az.ca.nus';
    var ons='';
    var adres='';
    for (i=name.length;i>=0;i--)
    {
      adres+=name.slice(i,i+1);
    }
    if (optdom!='')
    {
      site=optdom;
    }
    for (i=site.length;i>=0;i--)
    {
      ons+=site.slice(i,i+1);
    }
    adres += '@';
	var heading='<a onmouseover="window.status=\'' + adres + ons + '\'" onmouseout="window.status=document.title" HREF="';
	//var heading='<a onmouseover="window.status=\'' + htext + '\';return true" onmouseout="window.status=document.title" HREF="';
    output = heading + 'mailto:' + adres + ons + '">' + htext + '</a>';
    return output;
  }

function changeBkgColor(objName, colour) {
	//var colour1 = document.getElementById(objName).style.backgroundColor
	var obj = document.getElementById(objName)
    obj.style.backgroundColor = colour;
	document.window.status = objName
 }
