src=$1
grep "incl .*/"$src \
    $BASILISK/*.tags \
    $BASILISK/navier-stokes/*.tags \
    $BASILISK/examples/*.tags \
    $BASILISK/test/*.tags | \
    awk -v basilisk=$BASILISK '
            function title(fname) {
              getline <fname
              if ($1 != "/**")
                return " ";
              while ($1 != "#")
                getline <fname;
              gsub("# ", "", $0);
              gsub("*/$", "", $0);
              return $0;
            }
            {
              gsub(basilisk, "", $1);
              gsub(".tags:.*", "", $1);
              used = "/src" $1;
              lineno = $4;
              print "used " title(basilisk $1) "\t" used " " lineno;
            }'
