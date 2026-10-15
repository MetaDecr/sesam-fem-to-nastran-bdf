#!/usr/bin/env python3
"""fem2bdf GUI 服务（仅依赖 Python 标准库）。

启动后在浏览器打开页面，拖入 .fem 文件即可调用同目录上一级的
fem2bdf.exe 完成转换，返回 BDF 文本与转换日志。

    python server.py [--host 127.0.0.1] [--port 7100] [--open]

环境变量 FEM2BDF_EXE 可覆盖转换器路径。
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def _default_exe():
    """Locate fem2bdf.exe in dev layout and PyInstaller bundle alike."""
    if getattr(sys, "frozen", False):
        # PyInstaller one-file bundle: data/binaries live in _MEIPASS.
        return Path(sys._MEIPASS) / "fem2bdf.exe"
    return ROOT.parent / "fem2bdf.exe"


INDEX = (Path(sys._MEIPASS) / "index.html") if getattr(sys, "frozen", False) \
    else ROOT / "index.html"
EXE = Path(os.environ.get("FEM2BDF_EXE", str(_default_exe())))

MAX_UPLOAD = 512 * 1024 * 1024   # 512 MB
CONVERT_TIMEOUT = 600            # seconds


def safe_name(name, default="model.fem"):
    """Strip any path components and hostile characters from a file name."""
    name = os.path.basename(name or "").strip()
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    return name or default


class Handler(BaseHTTPRequestHandler):
    server_version = "fem2bdf-gui/1.0"

    # -- helpers ------------------------------------------------------
    def _send(self, code, body, ctype):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, code=200):
        self._send(code, json.dumps(obj).encode("utf-8"),
                   "application/json; charset=utf-8")

    def log_message(self, fmt, *args):  # keep console noise down
        sys.stderr.write("[gui] %s\n" % (fmt % args))

    # -- routes -------------------------------------------------------
    def do_GET(self):
        path = self.path.split("?", 1)[0]
        if path in ("/", "/index.html"):
            self._send(200, INDEX.read_bytes(), "text/html; charset=utf-8")
        elif path == "/api/health":
            self._json({"ok": True, "exe": str(EXE), "exe_exists": EXE.exists()})
        else:
            self._json({"ok": False, "error": "not found"}, 404)

    def do_POST(self):
        if self.path.split("?", 1)[0] != "/api/convert":
            return self._json({"ok": False, "error": "not found"}, 404)

        try:
            length = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            length = 0
        if length <= 0:
            return self._json({"ok": False, "error": "空文件"}, 400)
        if length > MAX_UPLOAD:
            return self._json({"ok": False, "error": "文件超过 512 MB 限制"}, 413)

        fname = safe_name(self.headers.get("X-Filename"))
        if not fname.lower().endswith(".fem"):
            return self._json({"ok": False, "error": "只接受 .fem 文件"}, 400)

        if not EXE.exists():
            return self._json(
                {"ok": False,
                 "error": "未找到转换器 fem2bdf.exe（%s），请先构建，见 README" % EXE},
                500)

        data = self.rfile.read(length)
        work = Path(tempfile.mkdtemp(prefix="fem2bdf_gui_"))
        try:
            inp = work / fname
            inp.write_bytes(data)
            outp = work / (inp.stem + ".bdf")
            try:
                proc = subprocess.run(
                    [str(EXE), str(inp), "-o", str(outp)],
                    capture_output=True, timeout=CONVERT_TIMEOUT)
            except subprocess.TimeoutExpired:
                return self._json({"ok": False, "error": "转换超时（10 分钟）"}, 504)

            log = (proc.stdout + proc.stderr).decode("utf-8", "replace")
            bdf = None
            if outp.exists():
                bdf = outp.read_text(encoding="utf-8", errors="replace")
            ok = proc.returncode == 0 and bdf is not None
            self._json({"ok": ok,
                        "returncode": proc.returncode,
                        "log": log,
                        "bdf": bdf,
                        "output_name": outp.name})
        finally:
            shutil.rmtree(work, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description="fem2bdf web GUI server")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=7100)
    ap.add_argument("--open", action="store_true", default=None,
                    help="启动后自动打开浏览器")
    ap.add_argument("--no-open", dest="open", action="store_false",
                    help="启动后不打开浏览器")
    args = ap.parse_args()
    # 打包版（双击启动）默认打开浏览器；脚本运行默认不打开。
    if args.open is None:
        args.open = bool(getattr(sys, "frozen", False))

    httpd = ThreadingHTTPServer((args.host, args.port), Handler)
    url = "http://%s:%d/" % (args.host, args.port)
    note = "" if EXE.exists() else "  [警告] 未找到 fem2bdf.exe"
    print("fem2bdf GUI 已启动: %s\n转换器: %s%s\n按 Ctrl+C 停止。" % (url, EXE, note))
    if args.open:
        try:
            webbrowser.open(url)
        except Exception:
            pass
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
