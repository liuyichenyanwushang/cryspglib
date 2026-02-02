#!/usr/bin/env python3
import re

def fix_3x3_arrays_in_file(filename):
    with open(filename, 'r') as f:
        content = f.read()

    # 模式匹配错误的3x3数组格式
    # static NAME: [[TYPE; 3]; 3] = [a, b, c],
    #     [d, e, f],
    #     [g, h, i],;

    lines = content.split('\n')
    output_lines = []
    i = 0

    while i < len(lines):
        line = lines[i]

        # 检查是否是3x3数组声明行
        if re.match(r'^\s*static\s+([A-Z_]+):\s*\[\[(i32|f64);\s*3\];\s*3\]\s*=\s*\[', line):
            # 提取数组名和类型
            match = re.match(r'^\s*static\s+([A-Z_]+):\s*\[\[(i32|f64);\s*3\];\s*3\]\s*=\s*\[', line)
            name = match.group(1)
            type_name = match.group(2)

            # 收集完整的数组定义（可能跨多行）
            array_def = line
            line_num = i + 1

            # 继续收集直到遇到分号
            while line_num < len(lines) and not lines[line_num-1].strip().endswith(';'):
                array_def += '\n' + lines[line_num]
                line_num += 1

            # 现在处理这个数组定义
            # 提取所有数字
            numbers = re.findall(r'[-]?\d+\.?\d*/?\d*\.?\d*', array_def)

            # 应该至少有9个数字
            if len(numbers) >= 9:
                # 格式化为正确的Rust 3x3数组
                row1 = f'[{numbers[0]}, {numbers[1]}, {numbers[2]}]'
                row2 = f'[{numbers[3]}, {numbers[4]}, {numbers[5]}]'
                row3 = f'[{numbers[6]}, {numbers[7]}, {numbers[8]}]'
                fixed_line = f'static {name}: [[{type_name}; 3]; 3] = [{row1}, {row2}, {row3}];'
                output_lines.append(fixed_line)

                # 跳过已处理的行
                i = line_num - 1
            else:
                # 无法修复，保留原样
                output_lines.append(line)
        else:
            output_lines.append(line)

        i += 1

    # 写回文件
    with open(filename, 'w') as f:
        f.write('\n'.join(output_lines))

    print(f"修复了文件: {filename}")

if __name__ == '__main__':
    fix_3x3_arrays_in_file('src/hall_symbol.rs')