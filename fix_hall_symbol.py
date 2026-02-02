#!/usr/bin/env python3
"""
修复hall_symbol.rs文件的语法错误和数组大小问题
"""

import re

def fix_3x3_arrays(content):
    """修复3x3数组的语法格式"""
    # 模式匹配：static NAME: [[TYPE; 3]; 3] = [a, b, c],
    # 应该改为：static NAME: [[TYPE; 3]; 3] = [[a, b, c], [d, e, f], [g, h, i]];

    lines = content.split('\n')
    fixed_lines = []
    i = 0

    while i < len(lines):
        line = lines[i]
        # 检查是否是3x3数组声明
        if re.match(r'^\s*static\s+[A-Z_]+:\s*\[\[(i32|f64);\s*3\];\s*3\]\s*=', line):
            # 收集数组元素行
            array_lines = [line]
            i += 1
            while i < len(lines) and not lines[i].strip().endswith(';'):
                array_lines.append(lines[i])
                i += 1
            if i < len(lines):
                array_lines.append(lines[i])  # 分号行

            # 解析数组元素
            array_text = '\n'.join(array_lines)
            # 提取数组名称和类型
            match = re.match(r'^\s*static\s+([A-Z_]+):\s*\[\[(i32|f64);\s*3\];\s*3\]\s*=\s*(.*)', array_text, re.DOTALL)
            if match:
                name, type_name, rest = match.groups()
                # 提取所有行中的数字
                numbers = []
                for l in array_lines[1:]:
                    # 匹配数字，包括负数和分数
                    nums = re.findall(r'[-]?\d+\.?\d*/?\d*\.?\d*', l)
                    numbers.extend(nums)

                # 应该有9个数字
                if len(numbers) >= 9:
                    # 格式化为3x3数组
                    row1 = f'[{numbers[0]}, {numbers[1]}, {numbers[2]}]'
                    row2 = f'[{numbers[3]}, {numbers[4]}, {numbers[5]}]'
                    row3 = f'[{numbers[6]}, {numbers[7]}, {numbers[8]}]'
                    fixed_line = f'static {name}: [[{type_name}; 3]; 3] = [{row1}, {row2}, {row3}];'
                    fixed_lines.append(fixed_line)
                else:
                    # 无法修复，保留原样
                    fixed_lines.extend(array_lines)
            else:
                fixed_lines.extend(array_lines)
        else:
            fixed_lines.append(line)
        i += 1

    return '\n'.join(fixed_lines)

def fix_large_array_sizes(content):
    """修复大型数组的大小声明"""
    # 从C文件读取实际大小
    with open('src_c/hall_symbol.txt', 'r') as f:
        c_content = f.read()

    # 查找C数组定义并计算大小
    array_sizes = {}

    # 模式匹配：static TYPE name[][3][9] = { ... };
    pattern = r'static\s+(?:int|double)\s+([a-zA-Z_]+)\[\]\[3\]\[9\]\s*=\s*\{([^{}]*(?:\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}[^{}]*)*)\};'

    # 更简单的模式：查找数组定义并计算最外层块
    lines = c_content.split('\n')
    current_array = None
    brace_count = 0
    in_array = False
    array_blocks = []

    for line in lines:
        # 查找数组定义开始
        if re.match(r'^\s*static\s+(?:int|double)\s+([a-zA-Z_]+)\[\]\[3\]\[9\]\s*=\s*\{', line):
            match = re.match(r'^\s*static\s+(?:int|double)\s+([a-zA-Z_]+)\[\]\[3\]\[9\]\s*=\s*\{', line)
            if match:
                current_array = match.group(1)
                brace_count = 1
                in_array = True
                array_blocks = []
                # 检查行中是否已经有更多大括号
                brace_count += line.count('{') - 1
                brace_count -= line.count('}')
        elif in_array:
            brace_count += line.count('{')
            brace_count -= line.count('}')

            # 当brace_count为1时，我们在最外层，可以计算块
            if brace_count == 1 and '}' in line:
                # 块结束
                pass

            if brace_count == 0:
                # 数组结束
                in_array = False
                # 计算大小（简化：计算顶层块数）
                # 对于这些数组，每个顶层块对应一个元素
                # 我们可以通过计算特定模式来估算
                pass

    # 已知的数组大小（从编译错误中收集）
    known_sizes = {
        'TRICLI_VSPU': 2,
        'TRICLI_GENERATORS': 2,
        'MONOCLI_VSPU': 13,
        'MONOCLI_GENERATORS': 13,
        'MONOCLI_A_VSPU': 8,
        'MONOCLI_B_VSPU': 8,
        'MONOCLI_C_VSPU': 8,
        'MONOCLI_I_VSPU': 8,
        'ORTHO_VSPU': 13,
        'ORTHO_GENERATORS': 13,
        'ORTHO_F_VSPU': 8,
        'ORTHO_I_VSPU': 8,
        'ORTHO_A_VSPU': 8,
        'ORTHO_B_VSPU': 8,
        'ORTHO_C_VSPU': 8,
        'TETRA_VSPU': 8,
        'TETRA_GENERATORS': 8,
        'TETRA_I_VSPU': 8,
        'TRIGO_VSPU': 13,
        'TRIGO_GENERATORS': 13,
        'RHOMBO_H_VSPU': 8,
        'RHOMBO_H_GENERATORS': 8,
        'RHOMBO_P_VSPU': 8,
        'RHOMBO_P_GENERATORS': 8,
        'HEXA_VSPU': 8,
        'HEXA_GENERATORS': 8,
        'CUBIC_VSPU': 13,
        'CUBIC_GENERATORS': 13,
        'CUBIC_F_VSPU': 8,
        'CUBIC_I_VSPU': 8,
    }

    # 修复Rust文件中的大小声明
    for array_name, size in known_sizes.items():
        # 查找模式：static ARRAY_NAME: [[[ TYPE; 9]; 3]; N] =
        pattern = rf'(static\s+{array_name}:\s*\[\[\[\s*(i32|f64);\s*9\];\s*3\];\s*)(\d+)(\]\s*=)'
        def repl(match):
            return f'{match.group(1)}{size}{match.group(4)}'

        content = re.sub(pattern, repl, content)

    return content

def main():
    # 读取生成的文件
    with open('src/hall_symbol.rs', 'r') as f:
        content = f.read()

    print("修复3x3数组格式...")
    content = fix_3x3_arrays(content)

    print("修复大型数组大小...")
    content = fix_large_array_sizes(content)

    # 写回文件
    with open('src/hall_symbol.rs', 'w') as f:
        f.write(content)

    print("修复完成！")

if __name__ == '__main__':
    main()