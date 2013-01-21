#include <stdio.h>
#include <string.h>

#define N 7

int main (int argc, char ** argv)
{
  if (argc != 2) {
    fprintf (stderr, "usage: endfor FILE\n");
    return 1;
  }
  FILE * fp = fopen(argv[1], "r");
  if (fp == NULL) {
    perror (argv[1]);
    return 1;
  }

  char c, buf[N+1] = "";
  char foreachs[80]; int nforeachs = 0;
  char endforeachs[80]; int nendforeachs = 0;
  int n = 0, comment = 0, comments = 0, para = 0, line = 1, foreach = 0, inforeach = 0;
  int endfor = 0, foreachlevel = 0;
  printf ("# 1 \"%s\"\n", argv[1]);
  while ((c = fgetc(fp)) != EOF) {
    putchar(c);
    if (foreach) {
      if ((c >= 'a' && c <= 'z') || c == '_')
	foreachs[nforeachs++] = c;
      else {
	foreachs[nforeachs++] = '\0';
	inforeach = 1; foreachlevel = para;
	foreach = 0;
      }
    }
    if (endfor) {
      if ((c >= 'a' && c <= 'z') || c == '_')
	endforeachs[nendforeachs++] = c;
      else {
	endforeachs[nendforeachs++] = '\0';
	if (strcmp(&endforeachs[4], foreachs)) {
	  fprintf (stderr, 
		   "%s:%d: error: "
		   "foreach%s() loop ended with end_for%s()\n", 
		   argv[1], line, foreachs, endforeachs);
	  return 1;
	}
	endfor = 0; nendforeachs = 0;
      }
    }

    switch (c) {
    case ';':
      if (!comment && !comments && inforeach && para == foreachlevel) {
	inforeach = 0;
	printf (" end_foreach%s(); ", foreachs);
      }
      break;
    case '*':
      if (buf[n-1] == '/')
	comment = 1;
      break;
    case '/':
      if (buf[n-1] == '*')
	comment = 0;
      else if (!comment && buf[n-1] == '/')
	comments = 1;
      break;
    case '\n': comments = 0; line++; break;
    case '{': 
      if (!comment && !comments) para++; 
      break;
    case '}': 
      if (!comment && !comments) {
	para--;
	if (para < 0) {
	  fprintf (stderr, 
		   "%s:%d: error: mismatched bracket\n", argv[1], line);
	  return 1;
	}
	if (inforeach && para == foreachlevel) {
	  inforeach = 0;
	  printf (" end_foreach%s(); ", foreachs);
	}
      }
      break;
    }
    if (n < N) {
      buf[n++] = c;
      buf[N] = '\0';
    }
    else {
      for (int i = 0; i < N - 1; i++)
	buf[i] = buf[i+1];
      buf[N-1] = c;
    }

    if (!comment && !comments) {
      if (!strncmp(buf, "end_for", 7)) {
	if (!inforeach) {
	  fprintf (stderr, "%s:%d: error: no loop to end\n", argv[1], line);
	  return 1;
	}
	inforeach = 0;
	endfor = 1;
      }
      else if (!endfor && !strncmp(buf, "foreach", 7)) {
	if (inforeach) {
	  fprintf (stderr, "%s:%d: error: foreach() loops cannot be nested\n", argv[1], line);
	  return 1;
	}
	foreach = 1; nforeachs = 0;
      }
    }
  }
  if (para != 0) {
    fprintf (stderr, "%s:%d: error: parse error at end of input\n", argv[1], line);
    return 1;
  }
}
