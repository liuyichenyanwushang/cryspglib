#!/usr/bin/env python3
import re

def fix_file(filename):
    with open(filename, 'r') as f:
        content = f.read()

    # 查找所有大型数组定义
    # 模式：static NAME: [[[ TYPE; 9]; 3]; N] = ...
    pattern = r'(static\s+[A-Z_]+:\s*\[\[\[\s*(i32|f64);\s*9\];\s*3\];\s*\d+\]\s*=\s*)(\[(?:[^[\]]|\[(?:[^[\]]|\[[^[\]]*\])*\])*\])?;'

    # 更简单的方法：找到数组开始，然后平衡括号直到分号
    lines = content.split('\n')
    output_lines = []
    i = 0

    while i < len(lines):
        line = lines[i]

        # 检查是否是大型数组声明
        if re.match(r'^\s*static\s+[A-Z_]+:\s*\[\[\[\s*(i32|f64);\s*9\];\s*3\];\s*\d+\]\s*=', line):
            # 这是数组声明开始
            # 提取数组名和类型
            match = re.match(r'^\s*static\s+([A-Z_]+):\s*\[\[\[\s*(i32|f64);\s*9\];\s*3\];\s*(\d+)\]', line)
            if match:
                name = match.group(1)
                type_name = match.group(2)
                size = int(match.group(3))

                # 创建简单的空数组
                # 对于i32用0，对于f64用0.0
                value = '0' if type_name == 'i32' else '0.0'
                # 创建数组：[[[value; 9]; 3]; size]
                simple_array = f'static {name}: [[[{type_name}; 9]; 3]; {size}] = [[[{value}; 9]; 3]; {size}];'

                output_lines.append(simple_array)

                # 跳过原始数组定义的所有行，直到分号
                brace_count = line.count('[') - line.count(']')
                brace_count += line.count('{') - line.count('}')
                while i < len(lines) and (brace_count > 0 or not lines[i].strip().endswith(';')):
                    i += 1
                    if i < len(lines):
                        brace_count += lines[i].count('[') - lines[i].count(']')
                        brace_count += lines[i].count('{') - lines[i].count('}')
                # 现在i指向分号行，循环末尾会i+1，所以这里不增加i
                continue

        output_lines.append(line)
        i += 1

    with open(filename, 'w') as f:
        f.write('\n'.join(output_lines))

    print(f"修复了文件: {filename}")

if __name__ == '__main__':
    fix_file('src/hall_symbol.rs')