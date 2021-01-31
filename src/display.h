/**
# Server-side display

This file implements the server-side of the interactive display of a
running Basilisk code. The client-side is typically done by the
[javascript bview](bview/README) implementation.

This works using two main components:

* [The wsServer WebSocket library](wsServer/README.md)
* [Vertex buffers](vertexbuffer.h)
*/

#ifndef DISPLAY_JS
# define DISPLAY_JS "http://basilisk.fr/three.js/editor/index.html"
#endif

#ifndef DISPLAY_HOST
# define DISPLAY_HOST "localhost"
#endif

#ifndef DISPLAY_RANGE
# define DISPLAY_RANGE "7100:7200"
#endif

#if 1
# define debug(...)
#else
# define debug(...) fprintf (stderr, __VA_ARGS__)
#endif

#include <wsServer/include/ws.h>
#pragma autolink -L$BASILISK/wsServer -lws

#include "view.h"
#include "khash.h"

typedef struct {
  int fd;
  int iter;
} DisplayClient;

KHASH_MAP_INIT_STR(strhash, DisplayClient *)

static struct {
  khash_t(strhash) * objects;
  int sock, port;
  char * error;
  Array * controls;
} Display = { .sock = -1 };

static void display_display()
{
  if (pid() == 0) {
    debug ("**************************\n");
    for (khiter_t k = kh_begin (Display.objects); k != kh_end (Display.objects);
	 ++k)
      if (kh_exist (Display.objects, k)) {
	debug ("%s", kh_key (Display.objects, k));
	DisplayClient * client = kh_value (Display.objects, k);
	while (client->fd >= 0) {
	  debug (" %d/%d", client->fd, client->iter);
	  client++;
	}
	debug ("\n");
      }
    debug ("--------------------------\n");
  }
}

static char * read_file_into_buffer (FILE * fp)
{
  if (fseek (fp, 0L, SEEK_END) < 0)
    return NULL;
  long bufsize = ftell (fp);
  if (bufsize <= 0)
    return NULL;
  char * buf = malloc (sizeof(char)*(bufsize + 1));
  if (fseek (fp, 0L, SEEK_SET) < 0) {
    free (buf);
    return NULL;
  }  
  size_t newLen = fread (buf, sizeof(char), bufsize, fp);
  buf[newLen] = '\0'; /* Just to be safe. */
  return buf;
}

static void display_command (const char * command)
{
  vertex_buffer_setup();
  VertexBuffer.vertex = 0; // fixme
  // Temporarily redirect stderr to catch errors
  fflush (stderr);
  int bak = dup (2);
  FILE * fp = tmpfile();
  dup2 (fileno (fp), 2);
  fclose (fp);

  char * line = strdup (command);
  process_line (line);
  free (line);
  
  free (Display.error);
  if (VertexBuffer.type >= 0)
    Display.error = NULL;
  else { // Nothing was drawn i.e. an error occured
    fflush (stderr);
    FILE * fp = fdopen (2, "r");
    Display.error = read_file_into_buffer (fp);
    int len = Display.error ? strlen (Display.error) : 0;
    if (len > 0 && Display.error[len - 1] == '\n')
      Display.error[len - 1] = '\0';
    fclose(fp);
  }
    
  // Restore stderr to its previous value
  dup2 (bak, 2);
  close (bak);

  if (VertexBuffer.type < 0)
    debug ("error: '%s'\n", Display.error);
  else
    debug ("position: %ld, normal: %ld, color: %ld, index: %ld, type: %d\n", 
	   VertexBuffer.position->len,
	   VertexBuffer.normal->len,
	   VertexBuffer.color->len,
	   VertexBuffer.index->len, VertexBuffer.type);
}

static int display_send (const char * command, int fd)
{
  debug ("sending '%s' to %d\n", command, fd);
  
  unsigned int commandlen = strlen (command);
  unsigned int errorlen = Display.error ? strlen (Display.error) : 0;

  int paddedlen = 4*ceil(commandlen/4.);
  size_t len = 2*sizeof(unsigned int) + paddedlen;

  if (errorlen > 0)
    len += errorlen;
  else
    len += 2*sizeof(int) +
      4*sizeof (unsigned int) +
      VertexBuffer.position->len +
      VertexBuffer.normal->len +
      VertexBuffer.color->len +
      VertexBuffer.index->len;
      
  if (ws_sendframe_init (fd, len, false, WS_FR_OP_BIN) < 0)
    return -1;
  
  if (ws_send (fd, &commandlen, sizeof(unsigned int)) < sizeof(unsigned int))
    return -1;
  if (ws_send (fd, &errorlen, sizeof(unsigned int)) < sizeof(unsigned int))
    return -1;
  if (ws_send (fd, command, commandlen) < commandlen)
    return -1;
  
  // padding to a multiple of four
  for (int i = 0; i < paddedlen - commandlen; i++) {
    char c = '\0';
    if (ws_send (fd, &c, 1) < 1)
      return -1;
  }

  if (errorlen > 0) {
    if (ws_send (fd, Display.error, errorlen) < errorlen)
      return -1;
  }
  else {
    if (ws_send (fd, &VertexBuffer.dim, sizeof(int)) < sizeof(int))
      return -1;
    if (ws_send (fd, &VertexBuffer.type, sizeof(int)) < sizeof(int))
      return -1;
    
    if (ws_send (fd, &VertexBuffer.position->len, sizeof (unsigned int))
	< sizeof (unsigned int))
      return -1;
    if (ws_send (fd, &VertexBuffer.normal->len, sizeof (unsigned int))
	< sizeof (unsigned int))
      return -1;
    if (ws_send (fd, &VertexBuffer.color->len, sizeof (unsigned int))
	< sizeof (unsigned int))
      return -1;
    if (ws_send (fd, &VertexBuffer.index->len, sizeof (unsigned int))
	< sizeof (unsigned int))
      return -1;
    
    if (ws_send (fd, VertexBuffer.position->p, VertexBuffer.position->len)
	< VertexBuffer.position->len)
      return -1;
    if (ws_send (fd, VertexBuffer.normal->p, VertexBuffer.normal->len)
	< VertexBuffer.normal->len)
      return -1;
    if (ws_send (fd, VertexBuffer.color->p, VertexBuffer.color->len)
	< VertexBuffer.color->len)
      return -1;
    if (ws_send (fd, VertexBuffer.index->p, VertexBuffer.index->len)
	< VertexBuffer.index->len)
      return -1;
  }
  
  return 0;
}

static void display_add (const char * command, int fd)
{
  debug ("adding '%s'\n", command);
  khiter_t k = kh_get (strhash, Display.objects, command);
  if (k == kh_end (Display.objects)) {
    int ret;
    k = kh_put (strhash, Display.objects, strdup (command), &ret);
    DisplayClient * client = malloc (sizeof(DisplayClient));
    client->fd = -1;
    kh_value (Display.objects, k) = client;
  }
  DisplayClient * clients = kh_value (Display.objects, k), * client = clients;
  int len = 0;
  while (client->fd >= 0) {
    if (client->fd == fd)
      fd = -1;
    client++, len++;
  }
  if (fd >= 0) { // not already in the array
    kh_value (Display.objects, k) = clients =
      realloc (clients, (len + 2)*sizeof (DisplayClient));
    clients[len].fd = fd;
    clients[len].iter = -1;
    clients[len + 1].fd = -1;
  }
  display_display();
}

#define JSON_BUILD(...) len += snprintf (build + len, 4096, __VA_ARGS__),	\
    build = realloc (build, len + 4096)

static void array_remove (khiter_t k, int fd)
{
  DisplayClient * clients = kh_value (Display.objects, k), * client = clients;
  int i = -1, len = 0;
  while (client->fd >= 0) {
    if (client->fd == fd) {
      if (i != -1)
	debug ("array_remove(): error! found multiple %d in '%s'\n",
	       fd, kh_key (Display.objects, k));
      i = len;
    }
    client++, len++;
  }
  if (i < 0)
    debug ("array_remove(): error! could not find %d in '%s'\n",
	   fd, kh_key (Display.objects, k));
  else if (len == 1) {
    free ((void *) kh_key (Display.objects, k));
    free ((void *) kh_value (Display.objects, k));
    kh_del (strhash, Display.objects, k);
  }
  else
    for (int j = i; j < len; j++)
      clients[j] = clients[j + 1];
}

static void display_remove (const char * command, int fd)
{
  debug ("removing '%s'\n", command);
  khiter_t k = kh_get (strhash, Display.objects, command);
  if (k == kh_end (Display.objects))
    debug ("display_remove(): error! could not find '%s' (%d)\n",
	   command, fd);
  else
    array_remove (k, fd);
  display_display();
}

typedef struct {
  char * name, * tooltip;
  void * ptr;
  double min, max;
  int size;
} DisplayControl;

static char * display_control_json()
{
  char * build = malloc (4096);
  int len = 0;
  JSON_BUILD ("#{");
  char sep = ' ';
  DisplayControl * d = Display.controls->p;
  for (int i = 0; i < Display.controls->len/sizeof(DisplayControl); i++, d++) {
    JSON_BUILD ("%c\n  \"%s\": { ", sep, d->name); sep = ',';
    if (d->tooltip)
      JSON_BUILD ("\"tooltip\": \"%s\", ", d->tooltip);
    switch (d->size) {
    case 4:
      JSON_BUILD ("\"type\": \"int\", \"value\": %d, \"min\": %g, \"max\": %g",
		  *((int *)d->ptr), d->min, d->max);
      break;
    case 8:
      JSON_BUILD ("\"type\": \"double\", \"value\": %g, \"min\": %g, \"max\": %g",
		  *((double *)d->ptr), d->min, d->max);
      break;
    default:
      assert (false);
    }
    JSON_BUILD (" }");
  }
  JSON_BUILD ("}");
  return build;
}

static DisplayControl * display_control_lookup (const char * name)
{
  DisplayControl * d = Display.controls->p;
  for (int i = 0; i < Display.controls->len/sizeof(DisplayControl); i++, d++)
      if (!strcmp (d->name, name))
	return d;
  return NULL;
}

struct _DisplayControl {
  void * ptr;
  double min, max;
  char * name, * tooltip, * ptr_name;
  int size;
};

static void display_control_internal (struct _DisplayControl p)
{
  if (Display.sock < 0) { // fixme: Display should be initialised in main()
    fprintf (stderr, "display_control(): error: cannot add control before initialisation\n");
    exit (1);
  }
  else {
    DisplayControl d;
    if (!p.name)
      p.name = p.ptr_name;

    if (display_control_lookup (p.name))
      return;
    
    d.name = strdup (p.name);
    d.tooltip = p.tooltip ? strdup (p.tooltip) : NULL;
    d.ptr = p.ptr;
    d.size = p.size;
    d.min = p.min, d.max = p.max;
    array_append (Display.controls, &d, sizeof (DisplayControl));
    
    char * controls = display_control_json();
    ws_sendframe_txt (0, controls, true);
    free (controls);
  }
}

#undef display_control
#define display_control(val, ...) display_control_internal (&(val), __VA_ARGS__, size = sizeof(val), ptr_name = #val)

static void display_control_update (const char * command, int fd)
{
  char * s = strdup (command), * s1 = strchr (command, ':');
  *s1++ = '\0';
  DisplayControl * d = display_control_lookup (command);
  if (d == NULL)
    debug ("display_control_update(): error! could not find '%s' (%d)\n",
	   command, fd);
  else {
    debug ("display_control_update (%s) = %s\n", command, s1);
    double val = atof(s1);
    if (d->max > d->min)
      val = clamp (val, d->min, d->max);
    switch (d->size) {
    case 4: *((int *)d->ptr) = val; break;
    case 8: *((double *)d->ptr) = val; break;
    default: assert (false);
    }

    char * controls = display_control_json();
    ws_sendframe_txt (- fd, controls, true);
    free (controls);
  }
  free (s);
}

static char * bview_interface_json()
{
  char * build = malloc (4096);
  int len = 0;

  JSON_BUILD ("{\n");
  
  char p[4096] = {0};
  int i = 0;
  while (bview_interface[i].json) {
    JSON_BUILD ("%s", i ? ",\n" : "");
    len += bview_interface[i].json (p, build + len, 4096);
    build = realloc (build, len + 4096);
    JSON_BUILD ("\n");
    i++;
  }
  JSON_BUILD ("}");
  return build;
}

void display_onclose (int fd)
{  
  debug ("closing %d\n", fd);
  for (khiter_t k = kh_begin (Display.objects); k != kh_end (Display.objects);
       ++k)
    if (kh_exist (Display.objects, k))
      array_remove (k, fd);
  display_display();
}

void display_onmessage (int fd, const char * msg, size_t size, int type)
{
  if (type == WS_FR_OP_TXT) {
    if (!msg)
      fprintf (stderr, "error receiving data on websocket\n");
    else switch (msg[0]) {
      case '+': display_add (msg + 1, fd); break;
      case '-': display_remove (msg + 1, fd); break;
      case '#': display_control_update (msg + 1, fd); break;
      default: fprintf (stderr,
			"display_onmessage: error: unknown message type '%s'\n",
			msg);
	break;
      }
  }
  else
    fprintf (stderr, "display_onmessage: error: unexpected message type '%d'\n",
	     type);
}

void display_onopen (int fd)
{
  char * interface = bview_interface_json();
  if (ws_sendframe_txt (fd, interface, false) < 0) {
    free (interface);
    display_onclose (fd);
    close (fd);
    return;
  }
  free (interface);
  
  char * controls = display_control_json();
  if (ws_sendframe_txt (fd, controls, false) < 0) {
    free (controls);
    display_onclose (fd);
    close (fd);
    return;
  }
  free (controls);

  if (display_defaults && ws_sendframe_txt (fd, display_defaults, false) < 0) {
    display_onclose (fd);
    close (fd);
  }
}

static void display_update (int i)
{
  for (khiter_t k = kh_begin (Display.objects); k != kh_end (Display.objects);
       ++k)
    if (kh_exist (Display.objects, k)) {
      DisplayClient * client = kh_value (Display.objects, k);
      while (client->fd >= 0) {
	if (client->iter < i)
	  break;
	client++;
      }
      if (client->fd >= 0) { // at least one client needs update
	const char * command = kh_key (Display.objects, k);
	display_command (command);
	client = kh_value (Display.objects, k);
	while (client->fd >= 0) {
	  if (client->iter < i) {
	      client->iter = i;
	      if (display_send (command, client->fd) < 0) {
		debug ("error sending '%s' to '%d'\n", command, client->fd);
		close (client->fd);
		display_onclose (client->fd);
		if (!kh_exist (Display.objects, k))
		  break;
	      }
	      else
		client++;
	  }
	  else
	    client++;
	}
	vertex_buffer_free();
	if (Display.error && kh_exist (Display.objects, k)) {
	  free ((void *) kh_key (Display.objects, k));
	  free ((void *) kh_value (Display.objects, k));
	  kh_del (strhash, Display.objects, k);
	}
      }
    }
}

int display_poll (int timeout)
{
  return ws_socket_poll (Display.sock,
			 display_onopen, display_onmessage, display_onclose,
			 timeout);
}

void display_url (FILE * fp)
{
  fprintf (fp, DISPLAY_JS "?ws://" DISPLAY_HOST ":%d", Display.port);
}

int display_usage = 20; // use 20% of runtime, maximum

/**
   -DDISPLAY=1, -DDISPLAY: play controls, start running immediately.
   -DDISPLAY=-1: play controls, initially paused.
   -DDISPLAY=0: no play controls, start running immediately.
  
   #include "display.h": play controls, initially paused.
   this can be changed by setting display_play = 0; in main().
*/

#ifdef DISPLAY
int display_play = DISPLAY < 0 ? -1 : 0; // negative: pause, zero: play, positive: step
#else
int display_play = -1;
#endif

static void display_destroy()
{
  if (pid() == 0) {
    for (khiter_t k = kh_begin (Display.objects); k != kh_end (Display.objects);
	 ++k)
      if (kh_exist (Display.objects, k)) {
	free ((void *) kh_key (Display.objects, k));
	free ((void *) kh_value (Display.objects, k));
      }
    kh_destroy (strhash, Display.objects);

    DisplayControl * d = Display.controls->p;
    for (int i = 0; i < Display.controls->len/sizeof(DisplayControl); i++, d++) {
      free (d->name);
      free (d->tooltip);
    }
    array_free (Display.controls);

    remove ("display.html");
  
    // fixme: close connection cleanly
    close (Display.sock);
  }
}

init_solver void display_init()
{
  if (pid() == 0) {
    const char * port = DISPLAY_RANGE;
    if (!strchr (port, ':'))
      Display.sock = ws_socket_open (atoi (port));
    else {
      char * s = strdup (port);
      char * s1 = strchr (s, ':');
      *s1++ = '\0';
      int pmax = atoi(s1);
      Display.port = atoi(s);
      while ((Display.sock = ws_socket_open (Display.port)) < 0 &&
	     Display.port < pmax)
	Display.port++;
      free (s);
    }
    if (Display.sock < 0) {
      char s[80];
      sprintf (s, "display(): could not open port '%s'", port);
      perror (s);
      exit (1);
    }
    Display.objects = kh_init (strhash);
    Display.controls = array_new();

    free_solver_func_add (display_destroy);

    FILE * fp = fopen ("display.html", "w");
    fputs ("<head><meta http-equiv=\"refresh\" content=\"0;URL=", fp);
    display_url (fp);
    fputs ("\"></head>\n", fp);
    fclose (fp);

#ifndef DISPLAY
    display_control (display_play, -1, 1, "Run/Pause");
#elif DISPLAY != 0
    display_control (display_play, -1, 1, "Run/Pause");
#endif
  
    display_control (display_usage, 0, 50, "Display %", 
		     "maximum % of runtime used by display");
  }
}

event refresh_display (i++, last)
{
  do {
    if (display_play)
      display_update (i);
    if (display_poll (display_play ? - 1 : 0))
      display_update (i);
  } while (display_play < 0);

  static timer global_timer = {0};
  static double poll_elapsed = 0.;
  if (poll_elapsed <= display_usage/100.*timer_elapsed (global_timer)) {
    global_timer = timer_start();
    display_update (i);
    poll_elapsed = timer_elapsed (global_timer);
  }
}
