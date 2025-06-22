import struct
import sys
with open(sys.argv[1], "rb") as f:
    data = f.read()
    # 假设你的结构体有一个 u32 和一个 f64，使用 struct.unpack
    values = struct.unpack("If", data)
    print(values)