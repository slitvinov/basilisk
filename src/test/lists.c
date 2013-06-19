// Check that lists work properly
// some compilers are picky about how lists of structures can be
// allocated

scalar sa[];
vector sb[];

scalar * list2 = {sa};
vector * list3 = {sb};

int main()
{
  scalar a = 1, b = 2, c = 3;
  vector d = {4, 5};
  scalar * list = {a,b,c,d};
  for (scalar s in list)
    fprintf (stderr, "%d\n", s);

  tensor e = {{6,7},{8,9}};
  vector * list1 = {d,e};
  for (vector v in list1)
    fprintf (stderr, "%d %d\n", v.x, v.y);
}
