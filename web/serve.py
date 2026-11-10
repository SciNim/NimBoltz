#!/usr/bin/env python3
"""Serves the page with cross-origin isolation, for threads:
    python3 web/serve.py [port]     # then open http://localhost:8000
Python's own `http.server` works too, with one thread."""
import functools, http.server, os, sys

class Handler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
        self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
        super().end_headers()

port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
handler = functools.partial(Handler, directory=os.path.dirname(os.path.abspath(__file__)))
http.server.ThreadingHTTPServer(('', port), handler).serve_forever()
