# -*- coding: utf-8 -*-
"""PPM(P6) -> PNG，纯标准库实现，供 GAMES101 作业看图用。

用法:
    python ppm2png.py <input.ppm> [output.png]

不带参数时默认读取同目录上层 out 构建目录里的 binary.ppm。
输出文件被图片查看器占用时会自动改名（result_时间戳.png），不会报错。
"""
import os
import sys
import time
import zlib
import struct


def read_ppm(path):
    f = open(path, "rb")
    buf = b""

    def nexttok():
        nonlocal buf
        while True:
            while buf and buf[:1] in b" \t\r\n":
                buf = buf[1:]
            if buf.startswith(b"#"):
                while buf and buf[:1] not in b"\r\n":
                    buf = buf[1:]
                continue
            if buf:
                i = 0
                while i < len(buf) and buf[i:i + 1] not in b" \t\r\n":
                    i += 1
                t = buf[:i]
                buf = buf[i:]
                return t
            b = f.read(65536)
            if not b:
                return None
            buf += b

    magic = nexttok()
    if magic != b"P6":
        raise ValueError("not a binary PPM (P6)")
    w = int(nexttok())
    h = int(nexttok())
    mx = int(nexttok())
    if mx != 255:
        raise ValueError("only 8-bit PPM supported (maxval=%d)" % mx)
    # maxval 之后恰好一个空白符要吃掉，否则整个像素数据错位 1 字节
    if buf.startswith(b"\r\n"):
        buf = buf[2:]
    elif buf[:1] in (b" ", b"\t", b"\r", b"\n"):
        buf = buf[1:]
    data = buf + f.read()
    f.close()
    if len(data) < w * h * 3:
        raise ValueError("truncated pixel data")
    return w, h, data[: w * h * 3]


def write_png(path, w, h, rgb):
    stride = w * 3
    raw = bytearray()
    for y in range(h):
        raw.append(0)  # filter type 0
        raw += rgb[y * stride:(y + 1) * stride]

    def chunk(t, d):
        c = t + d
        return struct.pack(">I", len(d)) + c + struct.pack(">I", zlib.crc32(c) & 0xFFFFFFFF)

    png = b"\x89PNG\r\n\x1a\n"
    png += chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
    png += chunk(b"IDAT", zlib.compress(bytes(raw), 6))
    png += chunk(b"IEND", b"")
    with open(path, "wb") as fo:
        fo.write(png)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "..", "out",
                                                             "build", "x64-Release",
                                                             "Assignment5", "binary.ppm")
    dst = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "result.png")

    if not os.path.exists(src):
        print("input not found:", src)
        return 1

    w, h, data = read_ppm(src)
    try:
        write_png(dst, w, h, data)
    except PermissionError:
        # 上一张还被图片查看器锁着，换个名字
        dst = os.path.splitext(dst)[0] + time.strftime("_%H%M%S") + ".png"
        write_png(dst, w, h, data)
    print("%dx%d -> %s" % (w, h, os.path.abspath(dst)))
    print("OUT=" + os.path.abspath(dst))
    return 0


if __name__ == "__main__":
    sys.exit(main())
