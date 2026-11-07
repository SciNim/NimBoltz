#!/usr/bin/env python3
"""Rewrites the single precision literals (`0.987'f32`) inside double
precision tables (`array[N, float64] = [...]`) of gas modules as the double
they widen to (0.9869999885559082). Nim folds `'f32` literals in a `const`
float64 array from the decimal text, without rounding to single precision
first; at run time (`let`) they are rounded. Usage: widen_f32.py file.nim..."""
import re, struct, sys

lit = re.compile(r"(-?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?)'f32")

def widened(m):
    v = struct.unpack('f', struct.pack('f', float(m.group(1))))[0]
    return repr(v)

for path in sys.argv[1:]:
    src = open(path).read()
    out, pos, n = [], 0, 0
    for m in re.finditer(r"(?:let|const) +\w+\*?: array\[[^\]]*float64\] = \[", src):
        end = src.index(']', m.end())
        body = src[m.end():end]
        new, k = lit.subn(widened, body)
        n += k
        if k:   # reflow to 80 columns
            items = [x.strip() for x in new.replace('\n', ' ').split(',') if x.strip()]
            lines, cur = [], '  '
            for it in items:
                if len(cur) + len(it) + 2 > 80 and cur.strip():
                    lines.append(cur.rstrip()); cur = '  '
                cur += it + ', '
            lines.append(cur.rstrip().rstrip(','))
            new = '\n' + '\n'.join(lines)
        out.append(src[pos:m.end()]); out.append(new); pos = end
    out.append(src[pos:])
    open(path, 'w').write(''.join(out))
    print(f"{path}: {n} literals widened")
