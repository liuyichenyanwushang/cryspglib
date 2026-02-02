#!/usr/bin/env python3
import re

with open('src/hall_symbol.rs', 'r') as f:
    lines = f.readlines()

# 找到 TRICLI_VSPU 的开始和 TRICLI_GENERATORS 的开始
start = -1
end = -1
for i, line in enumerate(lines):
    if 'static TRICLI_VSPU:' in line and start == -1:
        start = i
    if 'static TRICLI_GENERATORS:' in line and start != -1 and end == -1:
        end = i
        break

if start != -1 and end != -1:
    # 注释掉这些行
    for i in range(start, end):
        lines[i] = '// ' + lines[i]

    with open('src/hall_symbol.rs', 'w') as f:
        f.writelines(lines)
    print(f"注释了行 {start} 到 {end-1}")

    # 在注释后添加一个简单的数组定义
    simple_array = 'static TRICLI_VSPU: [[[f64; 9]; 3]; 2] = [[[0.0; 9]; 3]; 2];\n'
    lines.insert(start, simple_array)

    with open('src/hall_symbol.rs', 'w') as f:
        f.writelines(lines)
    print("添加了简单的数组定义")
else:
    print("找不到数组范围")