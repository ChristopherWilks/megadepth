f=$1
export LC_ALL=C
fgrep -v "		" $f | sort -u > $f.new
mv $f.new $f
