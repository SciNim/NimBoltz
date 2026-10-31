# Reference Fortran build and its cached outputs, shared by tests/regress
# and tests/stat. Source it, then:
#   reference_setup             builds $ref/magboltz-f if needed
#   reference_cmds deck...      prints one command per deck whose cached
#                               output $ref/<deck>.out is missing or stale
#                               ($ref/<deck>.rc: exit status, 124 = timeout)
# Env: MAGBOLTZ_SRC, WORK, FULL (1: unmodified NSCALE; default 4e7 -> 4e6),
#      TIMEOUT (seconds per run, default 3600).
src=${MAGBOLTZ_SRC:-$HOME/src/Magboltz/magboltz-11.21.f}
work=${WORK:-/tmp/nbregress}
to=${TIMEOUT:-3600}
mode=fast; [ "${FULL:-0}" = 1 ] && mode=full
nscale=4000000; [ $mode = full ] && nscale=40000000
ref=$work/ref-$mode
FFLAGS="-O0 -ffp-contract=off -std=legacy -finit-local-zero"

reference_setup() {
  mkdir -p "$ref"
  if [ $mode = full ]; then
    cp "$src" "$ref/magboltz.f.new"
  else
    sed 's/^      NSCALE=40000000 *$/      NSCALE=4000000/' "$src" > "$ref/magboltz.f.new"
    [ "$(grep -c '^      NSCALE=4000000$' "$ref/magboltz.f.new")" = 2 ] || { echo "NSCALE patch failed"; return 1; }
  fi
  cmp -s "$ref/magboltz.f.new" "$ref/magboltz.f" || mv "$ref/magboltz.f.new" "$ref/magboltz.f"
  rm -f "$ref/magboltz.f.new"
  if [ ! "$ref/magboltz-f" -nt "$ref/magboltz.f" ]; then
    echo "== gfortran"; gfortran $FFLAGS "$ref/magboltz.f" -o "$ref/magboltz-f"
  fi
}

reference_cmds() {
  local d n
  for d in "$@"; do
    n=$(basename "$d" .in)
    if [ ! "$ref/$n.out" -nt "$ref/magboltz-f" ] || [ ! "$ref/$n.out" -nt "$d" ]; then
      echo "timeout $to $ref/magboltz-f < $d > $ref/$n.out.tmp 2>/dev/null; echo \$? > $ref/$n.rc; mv $ref/$n.out.tmp $ref/$n.out"
    fi
  done
}
